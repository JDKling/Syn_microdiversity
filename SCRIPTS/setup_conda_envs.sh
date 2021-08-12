#!/usr/bin/env bash
test=('/CONDA-ENVS/env0.yml' '/CONDA-ENVS/env1.yml' '/CONDA-ENVS/env2.yml' '/CONDA-ENVS/env3.yml' '/CONDA-ENVS/env4.yml')

for yml in ${test[@]}; do
  if [ ! -f "$yml" ]; then
    echo "$yml  does not exist. Make sure all files got downloaded correctly."
    return 0
  fi
done
