#!/bin/bash
echo "=== pwd ==="; pwd
echo "=== whoami ==="; whoami
echo "=== HOME ==="; echo "$HOME"
echo "=== ls cwd ==="; ls -la
echo "=== build/ present? ==="; ls -la build 2>&1 | head
echo "=== ABE binary? ==="; ls -lh build/ABE build/TwoPunctureABE 2>&1
echo "=== /workspace ==="; ls -la /workspace 2>&1 | head
echo "=== is devpod home mounted? ==="; ls -la /home/h3250104945/HPC101/src/lab4/build 2>&1 | head
echo "=== nproc ==="; nproc
echo "=== DONE ==="
