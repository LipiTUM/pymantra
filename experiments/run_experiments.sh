#!/usr/bin/bash

# static variables
OK="\033[0;32m"
FAIL="\033[0;31m"
NC="\033[0m"

echo "====="
echo "starting mantra experiments"
echo -e "=====\n"

# metabolome-only enrichment
echo "Running metabolome analysis..."
if python3 pre_process_metabolome.py && python3 run_metabolome.py; then
  echo -e "${OK}Metabolome results computed${NC}"
else
   echo -e "${FAIL}Metabolome analysis failed${NC}"
fi

# metabolome-microbiome enrichment
echo -e "\nRunning metabolome-microbiome analysis..."
# if python3 pre_process_metabolome.py &&
#     python3 run_metabolome_microbiome_enrichment.py &&
#     python3 multiomics_downstream_analysis.py; then
#   echo -e "${OK}Metabolome-Microbiome results computed${NC}"
# else
#    echo -e "${FAIL}Metabolome-Microbiome analysis failed${NC}"
# fi

echo -e "\n====="
echo "mantra experiments finished"
echo "====="
