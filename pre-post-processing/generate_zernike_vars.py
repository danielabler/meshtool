import pandas as pd

data = pd.read_excel("/data/ownCloud/MY_PROJECTS/CorneaMeshing/new_cornea/CuencaEye.xlsx", skiprows=4) 
data.columns = ["ant_id", "ant_zernike","post_id", "post_zernike"]

with open("/data/ownCloud/MY_PROJECTS/CorneaMeshing/new_cornea/zernike.txt", 'w') as outfile:
    outfile.write("\n =========ANTERIOR SURFACE========= \n\n")    
    for index, row in data.iterrows():
        outfile.write("j_coeff_map_ant[%i] = %s;\n" %(int(row.ant_id)-1, row.ant_zernike))
        
    outfile.write("\n =========POSTERIOR SURFACE========= \n\n")    
    for index, row in data.iterrows():
        outfile.write("j_coeff_map_post[%i] = %s;\n" %(int(row.post_id)-1, row.post_zernike))    
        
        
with open("/data/ownCloud/MY_PROJECTS/CorneaMeshing/new_cornea/zernike_xml.txt", 'w') as outfile:
    outfile.write("\n =========ANTERIOR SURFACE========= \n\n")    
    for index, row in data.iterrows():
        outfile.write("<ZernikeCoefficientSingleIndex j='%d'>%.20f</ZernikeCoefficientSingleIndex>\n" %(int(row.ant_id)-1, row.ant_zernike))
        
    outfile.write("\n =========POSTERIOR SURFACE========= \n\n")    
    for index, row in data.iterrows():
        outfile.write("<ZernikeCoefficientSingleIndex j='%d'>%.20f</ZernikeCoefficientSingleIndex>\n" %(int(row.post_id)-1, row.post_zernike))
