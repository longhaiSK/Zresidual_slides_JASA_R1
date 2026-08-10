#!/bin/sh
set -e  # Stop immediately if any build fails

# Directories
SUBM_DIR="./submissions"
BLIND_DIR="./build/blinded"
UNBLIND_DIR="./build/unblinded"

mkdir -p "$SUBM_DIR" "$BLIND_DIR" "$UNBLIND_DIR"

# Cover letter: render with Quarto, then copy the PDF(s) to submissions
for f in covers/*.qmd; do
    quarto render "$f" --to typst
    cp "covers/$(basename "$f" .qmd).pdf" "$SUBM_DIR"
done

# build_pdf <source.tex> <build_dir> <anon flag: 0=blinded, 1=unblinded>
build_pdf() {
    src="$1"
    bdir="$2"
    anon="$3"

    echo "Building $src in $bdir (Anon=$anon)..."

    latexmk -pdf \
      -outdir="$bdir" \
      -interaction=nonstopmode \
      -usepretex -pretex="\def\anon{$anon}" \
      "$src"

    # Copy the newly generated .aux file to the root directory.
    # This allows \externaldocument to find it, and keeps editor autocomplete working.
    src_base=$(basename "$src" .tex)
    cp "$bdir/$src_base.aux" "./"
}

echo "=== BUILDING BLINDED VERSIONS ==="
# Pass 1: Main (Generates JASA-template.aux and copies it to root)
build_pdf JASA-template.tex "$BLIND_DIR" 0

# Pass 2: Supplement (Finds root JASA-template.aux, generates supplement.aux and copies to root)
build_pdf supplement.tex "$BLIND_DIR" 0

# Pass 3: Main (Finds root supplement.aux, resolves all references)
build_pdf JASA-template.tex "$BLIND_DIR" 0

# Rename and copy to submissions
cp "$BLIND_DIR/JASA-template.pdf" "$SUBM_DIR/zresdi_blinded.pdf"
cp "$BLIND_DIR/supplement.pdf" "$SUBM_DIR/supp_blinded.pdf"


echo "=== BUILDING UNBLINDED VERSIONS ==="
# Because the .aux files are still in the root, cross-references will 
# resolve instantly on the first pass here.
# Pass 1: Main
build_pdf JASA-template.tex "$UNBLIND_DIR" 1

# Pass 2: Supplement 
build_pdf supplement.tex "$UNBLIND_DIR" 1

# Pass 3: Main 
build_pdf JASA-template.tex "$UNBLIND_DIR" 1

# Rename and copy to submissions
cp "$UNBLIND_DIR/JASA-template.pdf" "$SUBM_DIR/zresdi_unblinded.pdf"
cp "$UNBLIND_DIR/supplement.pdf" "$SUBM_DIR/supp_unblinded.pdf"


echo "=== BUILDING REVIEWER RESPONSES ==="
for f in reviewer_*.tex; do
    latexmk -pdf -outdir="./build" "$f"
    cp "./build/$(basename "$f" .tex).pdf" "$SUBM_DIR"
done

echo "All PDFs successfully built and copied to $SUBM_DIR."