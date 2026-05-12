import pandas as pd
import xlsxwriter
import openpyxl
import os
from rdkit import Chem
from rdkit.Chem.Draw.MolDrawing import MolDrawing,DrawingOptions

from rdkit.Chem import FragmentCatalog
from rdkit.Chem import Draw
from PIL import Image

def drawAllInOneFigure():
    fName='FunctionalGroups_for_plot.txt'
    fparams = FragmentCatalog.FragCatParams(1,6,fName)
    fparams.GetNumFuncGroups()
    mols=[]
    for i in range(fparams.GetNumFuncGroups()):
        mols.append(fparams.GetFuncGroup(i))
    img=Draw.MolsToGridImage(mols,molsPerRow=16)
    
    img.save('output_image.png')

# 根据smiles编码生成图片，生成的图片挡在# FG_images文件夹中
def generation_images(data):
    draw = data.smiles.tolist()
    for i in draw:
        mol = Chem.MolFromSmiles(i)
        Draw.MolToFile(mol,f'./FG_images/img{i}.png',size=(150,100))

# 创建excel表格，将生成的图片读入
def load_images(data):
    workbook = xlsxwriter.Workbook('dataset_with_iamges.xlsx')
    worksheet = workbook.add_worksheet()
    for i,j in enumerate(data.smiles.tolist()):
        worksheet.write(f'A{i+1}', f'{j}')
        worksheet.insert_image(f'B{i+1}', f'./FG_images/img{j}.png')
    workbook.close()

def drawExcelWithImage():
    os.makedirs('./FG_images',exist_ok=True)
    # 读入要处理的数据
    df = pd.read_excel('dataset.xlsx', engine='openpyxl',header=None)
    generation_images(df)

def drawFGsInExcelWithImage():
    os.makedirs('./FG_images',exist_ok=True)
    fName='FunctionalGroups.txt'
    fparams = FragmentCatalog.FragCatParams(1,6,fName)
    fparams.GetNumFuncGroups()
    mols=[]
    names=[]
    for n in range(fparams.GetNumFuncGroups()):
        funcgroup= fparams.GetFuncGroup(n)
        mols.append(funcgroup)
        names.append(funcgroup.GetProp('_Name'))
    for m,mol in enumerate(mols):
        Draw.MolToFile(mol,f'./FG_images/img{m}.png',size=(150,100))
    workbook = xlsxwriter.Workbook('./dataset_with_images.xlsx')
    worksheet = workbook.add_worksheet()
    # 设置所有列的宽度
    worksheet.set_column(0, 10, 15)  # 设置第0列到第10列的宽度为15

    # 设置所有行的高度
    worksheet.set_default_row(80)  # 设置默认行高为30
    # worksheet.set_row(2, 60)  # 设置第2行的高度为60
    for i,j in enumerate(names):
        worksheet.write(f'A{i+1}', f'{j}')
        worksheet.insert_image(f'B{i+1}', f'./FG_images/img{i}.png', {'x_offset': 5, 'y_offset': 5, 'x_scale': 1, 'y_scale': 1})
    workbook.close()

if __name__ == '__main__':
    drawFGsInExcelWithImage()
    #drawAllInOneFigure()
