import xlrd

fname = 'method_AUC.xlsx'
fout = 'method_AUC.txt'
f2=open(fout,'w')

bk = xlrd.open_workbook(fname)
shxrange = range(bk.nsheets)
try:
	sh = bk.sheet_by_name("Sheet1")
except:
	print"no sheet"

nrows = sh.nrows
ncols = sh.ncols
print"nrows %d,ncols %d" %(nrows,ncols)
for i  in range(0,nrows):
	row_data = sh.row_values(i)
	f2.write(str(row_data[0]) +"&"+str(row_data[1]) +"&"+str(row_data[2]) +"&"+str(row_data[3]) +"&"+str(row_data[4]) +"&"+str(row_data[5]) +"&"+str(row_data[6]) +"\\"+"\\"+'\n')