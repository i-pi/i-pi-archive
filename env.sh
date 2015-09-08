ENV_BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export PATH=$ENV_BASE_DIR:$PATH/bin
export PYTHONPATH=$ENV_BASE_DIR:$PYTHONPATH

unset ENV_BASE_DIR
