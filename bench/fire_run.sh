#!/bin/sh
# FIRE's side of the reduction suite (bench/ibp.sh), run inside the pqcf-fire2 image (comparison/fire/Dockerfile) in
# /opt/fire/FIRE6, with the suite's directory at /out. For each line "Case Example Problem Target" of
# /out/fire_cases.txt: FIRE on its own example's start file, one thread, the target alone (timed:
# /out/<case>.fire.tables), then on PQCF's masters (/out/<case>.masters.m, untimed: /out/<case>.masters.tables).
# Times to /out/fire.json.
set -e
cd /opt/fire/FIRE6
config() {   # example, integrals file, output, config name
    V=$(sed -n 's/^#variables *//p' "examples/$1.config")
    P=$(sed -n 's/^#problem *//p' "examples/$1.config")
    PR=$(sed -n 's/^#preferred.*/&/p' "examples/$1.config")
    printf '#compressor none\n#threads 1\n#fthreads 1\n#calc flint\n#variables %s\n#start\n#folder examples/\n#problem %s\n#integrals %s\n%s\n#output %s\n' \
        "$V" "$P" "$2" "$PR" "$3" > "examples/$4.config"
}
echo "[" > /out/fire.json
SEP=""
while read -r CASE EX PR T; do
    printf '{{%s, {%s}}}\n' "$PR" "$T" > "examples/$CASE.m"
    config "$EX" "$CASE.m" "$CASE.tables" "$CASE"
    rm -rf temp
    S0=$(python3 -c 'import time; print(time.time())')
    bin/FIRE6 -c "examples/$CASE" > "/out/$CASE.fire.log" 2>&1
    SECS=$(python3 -c "import time; print(round(time.time() - $S0, 3))")
    cp "examples/$CASE.tables" "/out/$CASE.fire.tables"
    if [ -f "/out/$CASE.masters.m" ]; then
        cp "/out/$CASE.masters.m" "examples/$CASE.masters.m"
        config "$EX" "$CASE.masters.m" "$CASE.masters.tables" "$CASE.masters"
        rm -rf temp
        bin/FIRE6 -c "examples/$CASE.masters" > "/out/$CASE.masters.log" 2>&1
        cp "examples/$CASE.masters.tables" "/out/$CASE.masters.tables"
    fi
    printf '%s{"case": "%s", "fire_seconds": %s}\n' "$SEP" "$CASE" "$SECS" >> /out/fire.json
    SEP=","
    echo "FIRE $CASE: $SECS s"
done < /out/fire_cases.txt
echo "]" >> /out/fire.json
