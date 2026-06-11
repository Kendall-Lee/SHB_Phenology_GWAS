#!/bin/bash
# Run Enhanced QTL-Seq Analysis with BLUEs for all traits

cd /Users/kendalllee/Documents/Blueberry/LRLP_Blueberry/LRLP_QTLVar/DIY_QTLSeq

echo "============================================"
echo "Running QTL-Seq Analysis with BLUEs"
echo "============================================"
echo ""

# Check which phenotype files exist
echo "Checking for phenotype files..."
for trait in DTFruit Flow2Fruit FruitWT; do
  if [ -f "${trait}_SHB_allPheno.txt" ]; then
    echo "  ✓ Found: ${trait}_SHB_allPheno.txt"
  else
    echo "  ✗ Missing: ${trait}_SHB_allPheno.txt"
  fi
done

echo ""
echo "============================================"
echo ""

# Run analysis for each trait that has a phenotype file
for trait in DTFruit Flow2Fruit FruitWT; do
  if [ -f "${trait}_SHB_allPheno.txt" ]; then
    echo "Running analysis for: $trait"
    echo "--------------------------------------------"
    Rscript Enhanced_QTLSeq_Analysis_BLUES.R $trait

    if [ $? -eq 0 ]; then
      echo "✓ SUCCESS: $trait analysis complete"
    else
      echo "✗ ERROR: $trait analysis failed"
    fi
    echo ""
    echo "============================================"
    echo ""
  else
    echo "⊘ SKIPPING: $trait (phenotype file not found)"
    echo ""
  fi
done

echo "============================================"
echo "All analyses complete!"
echo "============================================"
echo ""
echo "Results saved in:"
for trait in DTFruit Flow2Fruit FruitWT; do
  if [ -d "results_${trait}" ]; then
    echo "  - results_${trait}/"
  fi
done
echo ""
