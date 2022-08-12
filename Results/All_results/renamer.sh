for file in  *.txt
do
	mv "$file" "${file/.txt/_WT.txt}"
done
