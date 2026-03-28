#!/usr/bin/env bash
# Short Shell Script to run Alphafold3 without writing the whole command everytime. Easily changeable for your purposes.
set -euo pipefail

SIF="/sybig/projects/tfclass/alphafold3.sif"
MODELS_DIR="/sybig/projects/tfclass/alphafold3_model_param"
DB_DIR="/sybig/projects/tfclass/public_databases"
CUDA_VISIBLE_DEVICES=0

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <input-json> <output-dir> [cuda-device]"
  echo "Example: $0 Your_Prediction.json Your_Prediction_OutputFolder 0"
  exit 1
fi

INPUT="$1"
OUTPUT="$2"
CUDA_VISIBLE_DEVICES="${3:-0}"

mkdir -p "$OUTPUT"

singularity exec --nv --env CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES" \
	--bind "$(dirname "$INPUT")":/root/af_input \
	--bind "$OUTPUT":/root/af_output \
	--bind "$MODELS_DIR":/root/models \
	--bind "$DB_DIR":/root/public_databases \
	"$SIF" python alphafold3/run_alphafold.py \
	--json_path=/root/af_input/"$(basename "$INPUT")" \
	--model_dir=/root/models \
	--db_dir=/root/public_databases \
	--output_dir=/root/af_output

