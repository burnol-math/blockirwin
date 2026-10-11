#!/usr/bin/env python3
"""Time and compare selected binary cases (default: 20/40/60 digits).

Run from the repository root. Each solver runs in its own subprocess;
subprocess timeouts are reported as unverified, never as agreement.
Recorded solver times exclude interpreter startup and imports; the
entire matrix may remain incomplete when the total budget expires.

    python benchmark_modalirwin_binary.py --digits 20 --seconds-per-case 8
    python benchmark_modalirwin_binary.py --digits 20,40,60 --total-seconds 170
"""
from __future__ import annotations
import argparse
import csv
import itertools
import json
from pathlib import Path
import subprocess
import sys
import time

WORDS = ('00', '11') + tuple(''.join(bits) for bits in itertools.product('01', repeat=3))
WORKER = '''import json, time
from {module} import {function}
start=time.perf_counter()
try:
    value={function}(2, {word!r}, {k}, {digits})
    print('RESULT_JSON:' + json.dumps({{'value':value,'seconds':time.perf_counter()-start}}),flush=True)
except Exception as exc:
    print('RESULT_JSON:' + json.dumps({{'exception':repr(exc),'seconds':time.perf_counter()-start}}),flush=True)
'''


def timed_call(module, function, word, k, digits, max_seconds):
    code = WORKER.format(module=module,function=function,word=word,k=k,digits=digits)
    try:
        done=subprocess.run([sys.executable,'-c',code], capture_output=True,
                            text=True, timeout=max_seconds, cwd=Path(__file__).parent)
    except subprocess.TimeoutExpired:
        return {'status': 'TIMEOUT', 'seconds': None}
    for line in done.stdout.splitlines():
        if line.startswith('RESULT_JSON:'):
            result=json.loads(line.removeprefix('RESULT_JSON:'))
            result['status']='OK' if 'value' in result else 'ERROR'
            return result
    return {'status':'ERROR', 'exception':done.stderr[-500:] or done.stdout[-500:]}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--digits',default='20,40,60',help='comma-separated precision values')
    ap.add_argument('--seconds-per-case',type=float,default=8,
                    help='maximum wall seconds for each single solver call')
    ap.add_argument('--total-seconds',type=float,default=170,
                    help='whole-script budget (below the 3-minute approval threshold)')
    ap.add_argument('--output',default='binary_modal_timings.csv')
    args=ap.parse_args()
    digits_values=[int(x) for x in args.digits.split(',')]
    deadline=time.monotonic()+args.total_seconds
    rows=[]
    output=Path(args.output)
    with output.open('w',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=(
            'base','word','k','digits','reference_status','reference_seconds',
            'reference_value','reference_error','modal_status','modal_seconds',
            'modal_value','modal_error','agree'))
        writer.writeheader()
        file.flush()
        for digits in digits_values:
            for word in WORDS:
                for k in (0,1):
                    if deadline-time.monotonic() <= 1:
                        print('Total time budget exhausted; unexecuted cases are not validated.',flush=True)
                        print(f'Wrote {len(rows)} attempted cases to {output}',flush=True)
                        return
                    row={'base':2,'word':word,'k':k,'digits':digits}
                    for module,function,label in [('blockirwin','blockirwin','reference'),
                                                  ('modalirwin','modalirwin','modal')]:
                        remaining=deadline-time.monotonic()
                        if remaining <= 1:
                            outcome={'status':'NOT_RUN','seconds':None}
                        else:
                            outcome=timed_call(module,function,word,k,digits,
                                              max(0.5,min(args.seconds_per_case,remaining-0.5)))
                        row[f'{label}_status']=outcome['status']
                        row[f'{label}_seconds']=outcome.get('seconds','')
                        row[f'{label}_value']=outcome.get('value','')
                        row[f'{label}_error']=outcome.get('exception','')
                    row['agree']=(row['reference_value']==row['modal_value']
                                  if row['reference_status']=='OK' and row['modal_status']=='OK'
                                  else 'UNVERIFIED')
                    rows.append(row)
                    writer.writerow(row)
                    file.flush()  # Preserve completed rows on interruption.
                    print(f"{word} k={k} digits={digits}: "
                          f"blockirwin={row['reference_status']}({row['reference_seconds']}); "
                          f"modalirwin={row['modal_status']}({row['modal_seconds']}); "
                          f"agree={row['agree']}",flush=True)
        print(f'Wrote {len(rows)} attempted cases to {output}',flush=True)


if __name__=='__main__':
    main()
