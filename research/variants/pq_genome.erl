%% Finding the variants in a sequenced genome, sliced (research/variants/README.md). The work is a list of numbered
%% tasks, each plain data, its result a pure function of it: alignment of one chunk of reads, then variant calling on
%% one slice of the genome. Tasks run on an executor (local slots, or one owner's Erlang nodes by erpc, the pattern
%% of pq_rexec), and the results are gathered in task order, so the output does not depend on how the work was
%% spread. The FVLLMP federation would be another executor for the same tasks.
-module(pq_genome).
-export([main/1, map/3, run_task/1, slices/3]).

%% {Name, Start, End} of the chromosome: numbered slices of Size bases; each task reads Margin bases more on either
%% side and keeps the variants at positions in [Start, End).
slices(Chrom, Len, Size) ->
    [{I, Chrom, S, min(S + Size, Len)} || {I, S} <- lists:zip(lists:seq(0, (Len - 1) div Size), lists:seq(0, Len - 1, Size))].

%%% Executors ---------------------------------------------------------------------------------------------------------

%% The results of Tasks in their order. Exec: local | {nodes, [node()]}; Slots: tasks at a time per node.
map(Tasks, local, Slots) ->
    map(Tasks, {nodes, [node()]}, Slots);
map(Tasks, {nodes, Nodes}, Slots) ->
    Indexed = lists:zip(lists:seq(1, length(Tasks)), Tasks),
    Free = lists:append(lists:duplicate(Slots, Nodes)),
    Done = loop(Indexed, Free, #{}, #{}, length(Tasks)),
    [maps:get(I, Done) || I <- lists:seq(1, length(Tasks))].

loop(_Queue, _Free, _Running, Done, N) when map_size(Done) =:= N -> Done;
loop([{I, T} | Queue], [Node | Free], Running, Done, N) ->
    Self = self(),
    Pid = spawn(fun() -> Self ! {done, self(), call(Node, T)} end),
    loop(Queue, Free, Running#{Pid => {I, T, Node}}, Done, N);
loop(Queue, Free, Running, Done, N) ->
    receive
        {done, Pid, {node_failed, Reason}} ->
            {{I, T, Node}, Running1} = maps:take(Pid, Running),
            io:format("task ~p: node ~p failed: ~p; retried~n", [I, Node, Reason]),
            loop([{I, T} | Queue], Free ++ [Node], Running1, Done, N);
        {done, Pid, {task_failed, Reason}} ->
            {{I, _T, Node}, _} = maps:take(Pid, Running),
            error({task_failed, I, Node, Reason});
        {done, Pid, {ok, Result}} ->
            {{I, _T, Node}, Running1} = maps:take(Pid, Running),
            loop(Queue, [Node | Free], Running1, Done#{I => Result}, N)
    end.

%% A task's own failure is the same on any node, so it ends the run (as in pq_rexec); only a node's failure is retried.
call(Node, T) ->
    try {ok, erpc:call(Node, ?MODULE, run_task, [T], infinity)}
    catch error:{exception, Reason, _Stack} -> {task_failed, Reason};
          error:{erpc, Reason} -> {node_failed, Reason}
    end.

%%% The tasks -----------------------------------------------------------------------------------------------------------

sh(Cmd) ->
    T0 = erlang:monotonic_time(millisecond),
    Port = open_port({spawn, "sh -c '" ++ Cmd ++ "'"}, [exit_status, stderr_to_stdout, binary]),
    Out = collect(Port, <<>>),
    {Out, erlang:monotonic_time(millisecond) - T0}.

collect(Port, Acc) ->
    receive
        {Port, {data, D}} -> collect(Port, <<Acc/binary, D/binary>>);
        {Port, {exit_status, 0}} -> Acc;
        {Port, {exit_status, S}} -> error({command_failed, S, Acc})
    end.

%% one chunk of read pairs: aligned with minimap2 (read group of the sample), sorted, indexed
run_task({align, Id, R1, R2, Ref, Out, Threads, Sample}) ->
    {_, Ms} = sh(io_lib:format("minimap2 -ax sr -t ~b -R \"@RG\\tID:chunk~b\\tSM:~s\" ~s ~s ~s 2>/dev/null"
                               " | samtools sort -@ 4 -m 1G -o ~s - 2>/dev/null && samtools index ~s",
                               [Threads, Id, Sample, Ref, R1, R2, Out, Out])),
    {align, Id, Out, Ms};
%% one slice: bcftools mpileup over every chunk's alignments in the slice and its margins, call, normalise, then keep
%% the variants whose position lies in the slice's own range. norm -c w: GRCh38 has a few IUPAC codes in its sequence
%% (e.g. B at chr3:16902879), where norm would otherwise stop; there it warns and leaves the record as called.
run_task({call, Id, Chrom, Start, End, Margin, Bams, Ref, Out}) ->
    R = io_lib:format("~s:~b-~b", [Chrom, max(1, Start + 1 - Margin), End + Margin]),
    {_, Ms} = sh(io_lib:format("bcftools mpileup -Ou -f ~s -r ~s -a AD,DP ~s 2>/dev/null"
                               " | bcftools call -mv -Ou 2>/dev/null | bcftools norm -c w -f ~s -Ou 2>/dev/null"
                               " | bcftools view -Oz -o ~s -t ~s:~b-~b 2>/dev/null && tabix -f -p vcf ~s",
                               [Ref, R, lists:join(" ", Bams), Ref, Out, Chrom, Start + 1, End, Out])),
    {call, Id, Out, Ms}.

%%% A run: chunks of reads in, one VCF out -----------------------------------------------------------------------------

%% main(Spec): Spec is the path of a JSON file:
%%   {"exec": "local" | ["node@host", ...], "align_slots": 1, "call_slots": 8, "threads": 16,
%%    "ref": ".../ref.fna", "index": ".../ref.mmi", "sample": "HG002", "chunks": [[r1, r2], ...],
%%    "regions": [["chr1", 248956422], ...] (or "chrom" and "length" for one), "slice": 1000000, "margin": 1000,
%%    "work": "/work", "out": "calls.vcf.gz", "skip_align": false}
%% Slices are numbered across the regions in their order, so the gathered calls follow the regions' order.
main(SpecPath) ->
    S = json:decode(element(2, file:read_file(SpecPath))),
    Exec = case maps:get(<<"exec">>, S) of <<"local">> -> local; Ns -> {nodes, [binary_to_atom(N) || N <- Ns]} end,
    [connect(N) || {nodes, Ns} <- [Exec], N <- Ns],
    W = binary_to_list(maps:get(<<"work">>, S)),
    Ref = binary_to_list(maps:get(<<"ref">>, S)),
    Chunks = maps:get(<<"chunks">>, S),
    Bams = [filename:join(W, io_lib:format("chunk~3..0b.bam", [I])) || I <- lists:seq(0, length(Chunks) - 1)],
    T0 = erlang:monotonic_time(millisecond),
    Aligned = case maps:get(<<"skip_align">>, S, false) of
        true -> [];
        false ->
            AlignTasks = [{align, I, binary_to_list(R1), binary_to_list(R2), binary_to_list(maps:get(<<"index">>, S)),
                           B, maps:get(<<"threads">>, S), binary_to_list(maps:get(<<"sample">>, S))}
                          || {I, [R1, R2], B} <- lists:zip3(lists:seq(0, length(Chunks) - 1), Chunks, Bams)],
            map(AlignTasks, Exec, maps:get(<<"align_slots">>, S))
    end,
    T1 = erlang:monotonic_time(millisecond),
    Regions = case maps:find(<<"regions">>, S) of
                  {ok, Rs} -> [{binary_to_list(C), L} || [C, L] <- Rs];
                  error -> [{binary_to_list(maps:get(<<"chrom">>, S)), maps:get(<<"length">>, S)}]
              end,
    Size = maps:get(<<"slice">>, S),
    Slices0 = lists:append([[{C, St, En} || {_, _, St, En} <- slices(C, L, Size)] || {C, L} <- Regions]),
    Slices = [{I, C, St, En} || {I, {C, St, En}} <- lists:zip(lists:seq(0, length(Slices0) - 1), Slices0)],
    CallTasks = [{call, I, C, St, En, maps:get(<<"margin">>, S), Bams, Ref,
                  filename:join(W, io_lib:format("slice~5..0b.vcf.gz", [I]))} || {I, C, St, En} <- Slices],
    Called = map(CallTasks, Exec, maps:get(<<"call_slots">>, S)),
    T2 = erlang:monotonic_time(millisecond),
    Out = binary_to_list(maps:get(<<"out">>, S)),
    List = Out ++ ".slices",
    ok = file:write_file(List, [[F, $\n] || {call, _, F, _} <- Called]),
    {_, _} = sh(io_lib:format("bcftools concat -Oz -o ~s -f ~s 2>/dev/null && tabix -f -p vcf ~s", [Out, List, Out])),
    T3 = erlang:monotonic_time(millisecond),
    Report = #{<<"align_seconds">> => (T1 - T0) / 1000, <<"call_seconds">> => (T2 - T1) / 1000,
               <<"gather_seconds">> => (T3 - T2) / 1000, <<"chunks">> => length(Chunks), <<"slices">> => length(Slices),
               <<"align_task_seconds">> => [Ms / 1000 || {align, _, _, Ms} <- Aligned],
               <<"call_task_seconds">> => [Ms / 1000 || {call, _, _, Ms} <- Called]},
    ok = file:write_file(Out ++ ".timing.json", json:encode(Report)),
    io:format("~s~n", [json:encode(Report)]).

connect(Node) ->
    true = net_kernel:connect_node(Node),
    {module, ?MODULE} = erpc:call(Node, code, ensure_loaded, [?MODULE]).
