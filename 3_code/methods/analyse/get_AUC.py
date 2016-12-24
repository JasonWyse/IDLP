import xlrd

fname = 'method_AUC.xlsx'
fout = 'method_AUC.txt'
f2=open(fout,'w')

bk = xlrd.open_workbook(fname)
shxrange = range(bk.nsheets)
try:
	sh = bk.sheet_by_name("Sheet2")
except:
	print"no sheet"

nrows = sh.nrows
ncols = sh.ncols
print"nrows %d,ncols %d" %(nrows,ncols)
for i  in range(0,nrows):
	row_data = sh.row_values(i)
	if(row_data[1]):
		f2.write( "%s&%.4f&%0.4f&%0.4f&%0.4f&%0.4f&%0.4f&%0.4f\\\\\n" % (row_data[0],row_data[1],row_data[2],row_data[3],row_data[4],row_data[5],row_data[6],row_data[7]))
	else:
		f2.write("%s&-&-&-&-&-&-&-\\\\\n" %(row_data[0]))
