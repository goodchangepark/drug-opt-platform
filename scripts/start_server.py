import subprocess
import time
import os
import sys
import urllib.request

def start():
    subprocess.run("fuser -k 8765/tcp", shell=True, capture_output=True)
    time.sleep(1)

    uvicorn_bin = "/home/xavier/chem/drug-opt-platform/.venv/bin/uvicorn"
    log_file = open("/home/xavier/chem/drug-opt-platform/uvicorn.log", "w")
    proc = subprocess.Popen(
        [uvicorn_bin, "backend.main:app", "--host", "127.0.0.1", "--port", "8765"],
        stdout=log_file,
        stderr=log_file,
        start_new_session=True,
        cwd="/home/xavier/chem/drug-opt-platform"
    )
    print(f"Spawned uvicorn with PID {proc.pid}")

    for i in range(25):
        time.sleep(1)
        try:
            with urllib.request.urlopen("http://127.0.0.1:8765/api/help/registry", timeout=2) as resp:
                if resp.status == 200:
                    print(f"Server is LIVE and healthy on port 8765 (waited {i+1}s)!")
                    return True
        except Exception:
            pass

    print("Failed to start uvicorn within 25s. Log contents:")
    log_file.close()
    with open("/home/xavier/chem/drug-opt-platform/uvicorn.log") as f:
        print(f.read())
    return False

if __name__ == "__main__":
    start()
