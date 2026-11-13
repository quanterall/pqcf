%% A pq_task for the tests (pq_task_tests): squares a number, declares its resources and describes itself.
-module(pq_task_example).
-behaviour(pq_task).
-export([run/2, resources/1, describe/1]).

run({square, N}, #{home := _}) -> N * N;
run({fail, Why}, _) -> error(Why).

resources({square, _}) -> #{cpu => 1};
resources(_) -> #{}.

describe({Kind, _}) -> #{kind => Kind}.
