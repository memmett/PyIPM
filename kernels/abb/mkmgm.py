
from xlwt import Workbook
from utils.io import read_csv

from glob import glob


def create_mgm_stand_sheet(sheet, name, standage, year, region=1, subregion=3):
    sheet.write(0, 0, 'Stand')
    sheet.write(0, 1, 'Stand ID')
    sheet.write(0, 2, name)

    sheet.write(1, 1, 'DateRun')
    sheet.write(1, 2, 'DUNNO')

    sheet.write(2, 1, 'StandAge')
    sheet.write(2, 2, standage)

    sheet.write(3, 1, 'Year')
    sheet.write(3, 2, year)

    sheet.write(4, 1, 'Region')
    sheet.write(4, 2, region)

    sheet.write(5, 1, 'SubRegion')
    sheet.write(5, 2, subregion)

    sheet.write(6, 1, 'CropPlanID')
    sheet.write(7, 1, 'CropPlansID')

    sheet.write(9,  2, 'Site')
    sheet.write(10, 1, 'All')
    sheet.write(11, 1, 'Con')
    sheet.write(12, 1, 'Dec')
    sheet.write(13, 1, 'Sw')
    sheet.write(13, 2, 18)
    sheet.write(14, 1, 'Pl')
    sheet.write(14, 2, 'NA')
    sheet.write(15, 1, 'Aw')
    sheet.write(15, 2, 20)
    sheet.write(16, 1, 'Sb')
    sheet.write(16, 2, 'NA')

    sheet.write(19, 0, 'Species')
    sheet.write(19, 1, 'Diameter (cm)')
    sheet.write(19, 2, 'Trees/ha')
    sheet.write(19, 3, 'Height (m)')
    sheet.write(19, 4, 'Total Age')
    sheet.write(19, 5, 'BHAge')


###############################################################################

book = Workbook()


for plotfile in glob('*.csv'):
    plot, _ = plotfile.split('.')
    meas    = read_csv(plotfile)

    sheet = book.add_sheet(plot)
    mgm_sheet(sheet, plot, 0, 0)

    for i, r in enumerate(meas):
        sheet.write(20+i, 0, r.species)         # species
        sheet.write(20+i, 1, float(r.dbh)/10.0) # diameter
        sheet.write(20+i, 2, 10.0)              # trees/ha
        sheet.write(20+i, 3, 10.0)              # height
        sheet.write(20+i, 4, 10.0)              # total age
        sheet.write(20+i, 5, 10.0)              # bhage


book.save('MGM Stands.xls')