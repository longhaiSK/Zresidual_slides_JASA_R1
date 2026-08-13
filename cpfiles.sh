# Set source and destination variables for cleaner commands
SRC="$Github/Zresid_Bayesian_JASA_R1"
DEST="."

# 1. Image and Plot Directories (from your original list)
cp -r $SRC/plot/ $DEST
cp -r $SRC/slidesplot $DEST

# 2. Slide Source Files (.qmd, .html, and associated _files folders)
cp -a $SRC/slides_* $DEST

# 3. Quarto Project Configuration (if you use project-level settings)
cp $SRC/*.yml $DEST


# 5. Bibliography and Styling (Modify if your extensions differ)
# Copies your .bib references and Reveal.js custom styling if they exist
cp $SRC/*.bib $DEST 2>/dev/null
cp $SRC/*.scss $DEST 2>/dev/null
cp $SRC/*.css $DEST 2>/dev/null

# 6. Saved R Data (Required if your slides load pre-computed Z-residual models)
cp -r $SRC/data $DEST 2>/dev/null
cp $SRC/*.rds $DEST 2>/dev/null
cp $SRC/*.RData $DEST 2>/dev/null
