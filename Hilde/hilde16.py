#
# This file is based on the file based on apcw_fracture.py
# 2016-03-15
#

### IMPORT PACKAGES AND DEFAULTS ###
print "Importing packages..."
import os
import sys
import re
sys.path.append("C:\\Program Files\\BEAMER\\v5.1.1_x64\\")
import BEAMERpy as beamer       # this is the official beamer python package.
from BEAMERpyDefaults import *  # this is my defaults file.

### INSTANTIATE ###
print "Instantiating BEAMER..."
lb = beamer.GBEAMER()
lb.set_number_of_threads(8)
lb.set_temp_file_behaviour({'Mode' : 'FixedDirectory',
                            'Directory' : os.getcwd()})

### SET DEFAULTS ###
print "Setting defaults..."
lb.extract_defaults(extract_def)
lb.bias_defaults(bias_def)
lb.fracture_defaults(fracture_def)
lb.fda_defaults(fda_def)
lb.transform_defaults(transform_def)
lb.export_gpf_defaults(export_gpf_def)
lb.import_gpf_defaults(import_gpf_def)

### IMPORT DXF ###
print "Importing DXF file..."
fn_in = 'hilde16.dxf'       # name of the dxf file to import
fn_prefix = re.sub(r'.dxf',r'',fn_in)
rotate_bln = 'top'

# in the future, this should probably just look in the folder for the only .dxf file.
dxf = lb.import_dxf({'LayerSet' : '*',
                     'DXFUnits' : 'um',
                     'DXFPolyMode' : 'ConvertToPolygon',
                     'MaxErrorForConversion' : str(import_res),
                     'SnappingRange' : str(import_res),
                     'FileName' : fn_in } )
## ROTATE ##
if rotate_bln == 'top':
    dxf = lb.transform(dxf,{'Rotation' : 90})
elif rotate_bln == 'bottom':
    dxf = lb.transform(dxf,{'Rotation' : -90})

### EXTRACT LAYERS ###
print "Extracting layers..."
cover = lb.extract(dxf,{'LayerSet' : 'cover'})
grid = lb.extract(dxf,{'LayerSet' : 'grid'})
holes = lb.extract(dxf,{'LayerSet' : 'holes'})
markers = lb.extract(dxf,{'LayerSet' : 'marks,markers'})
nitride = lb.extract(dxf,{'LayerSet' : 'nitride'})
phc = lb.extract(dxf,{'LayerSet' : 'phc'})
negY = lb.extract(dxf,{'LayerSet' : 'negYjunc'})
phccover = lb.extract(dxf,{'LayerSet' : 'phccover'})
prebias = lb.extract(dxf,{'LayerSet' : 'prebias'})
regionLayer = lb.extract(dxf,{'LayerSet' : 'regionLayer1,regionLayer2,regionLayer3,regionLayer4,regionLayer5,regionLayer6,regionLayer7,regionLayer8,regionLayer9,regionLayer10,regionLayer11,regionLayer12,regionLayer13,regionLayer14,regionLayer15,regionLayer16'})
regionLayer = lb.heal(regionLayer,{'TargetLayer' : 'regionLayer'})
tethers = lb.extract(dxf,{'LayerSet' : 'fda1,fda2,fda3,fda4,fda5,fda6,fda7,fda8'})
indicators = lb.extract(dxf,{'LayerSet' : 'indicators'})
vgroove = lb.extract(dxf,{'LayerSet' : 'vgroove'})
vgrooveFDA = lb.extract(dxf,{'LayerSet' : 'vgrooveFDA'})
window = lb.extract(dxf,{'LayerSet' : 'window'})

### BIASING ###
def coarseExp():
    '''Coarse exposure.'''
    windowMinusIndicators = lb.minus(window,indicators)
    windowMinusGrid = lb.minus(windowMinusIndicators,grid)
    vgrooveMinusGrid = lb.minus(vgroove,grid)
    biasPlus2prebias = lb.bias(lb.bool_or(prebias,cover),{'Bias' : 2})
    biasMinus2window = lb.bias(windowMinusIndicators,{'Bias' : -2})
    windowOutline = lb.minus(windowMinusIndicators,biasMinus2window)
    out = lb.merge([windowMinusGrid,vgrooveMinusGrid]) 
    # out is window and vee groove minus grid...
    # out = lb.merge([out,windowOutline,vgrooveFDA,markers])
    out = lb.merge([out,windowOutline,vgrooveFDA])
    # ... and the window outline, vee groove FDA and markers
    out = lb.minus(out,biasPlus2prebias)
    # ... but not nitride ...
    return out

def fineExp():
    '''Fine exposure.'''
    biasedPrebias = lb.bias(prebias, {'Bias' : 3})
    prebiasOutline = lb.minus(biasedPrebias,nitride) # exp around nitride
    biasedVgroove = lb.bias(vgroove, {'Bias' : -2})
    vgrooveOutline = lb.minus(vgroove,biasedVgroove) # exp around vgroove
    vgrooveOutline = lb.bool_and(vgrooveOutline,lb.bool_not(prebias))
    out = lb.bool_or(lb.minus(prebiasOutline,lb.bias(cover,{'Bias' : 3})),vgrooveOutline)
    out = lb.bool_and(out,vgroove) 
    # CHANGE THIS LATER TO FRACTURE
    # EACH ACCORDING TO ITS SHAPE
    return out

print "Compiling coarse exposure..."
coarse = coarseExp()

print "Compiling fine exposure..."
fine = fineExp()
phcBias = lb.merge([lb.minus(lb.minus(lb.bias(cover,{'Bias' : 3}),cover),nitride)])
fine = lb.merge([fine,phcBias])

# phc = lb.merge([phc,lb.bool_and(lb.minus(lb.bias(cover,{'Bias' : 3}),nitride),lb.bool_not(phccover))])
# phc = lb.merge([phc,negY])

### INTERMEDIATE EXPORT ###
print "Fracturing and performing intermediate export..."
lb.export_gpf(coarse, {'FileName' : fn_prefix+'_coarse.gpf',
                       'GridResolution' : coarse_res,
                       'SubfieldSizeX' : coarse_subfield_size,
                       'SubfieldSizeY' : coarse_subfield_size,
                       'MainfieldSizeX' : coarse_mainfield_size,
                       'MainfieldSizeY' : coarse_mainfield_size,
                       'BeamStepSize' : coarse_bss,
                       'MainFieldResolution' : coarse_mainfield_res,
                       'SubFieldResolution' : coarse_subfield_res,
                       'FractureMode' : 'LRFT'})

lb.export_gpf(lb.merge([fine,regionLayer]), 
              {'FileName' : fn_prefix+'_fine.gpf',
               'GridResolution' : fine_res,
               'SubfieldSizeX' : fine_subfield_size,
               'SubfieldSizeY' : fine_subfield_size,
               'MainfieldSizeX' : fine_mainfield_size,
               'MainfieldSizeY' : fine_mainfield_size,
               'BeamStepSize' : fine_bss,
               'MainFieldResolution' : fine_mainfield_res,
               'SubFieldResolution' : fine_subfield_res,
               'FractureMode' : 'Conventional',
               'RegionLayerSet' : 'regionLayer',
               'MainfieldPlacement' : 'RegionLayer'})

lb.export_gpf(lb.merge([tethers,regionLayer]),
              {'FileName' : fn_prefix+'_tethers.gpf',
               'GridResolution' : fine_res,
               'SubfieldSizeX' : fine_subfield_size,
               'SubfieldSizeY' : fine_subfield_size,
               'MainfieldSizeX' : fine_mainfield_size,
               'MainfieldSizeY' : fine_mainfield_size,
               'BeamStepSize' : fine_bss,
               'MainFieldResolution' : fine_mainfield_res,
               'SubFieldResolution' : fine_subfield_res,
               'FractureMode' : 'LRFT',
               'RegionLayerSet' : 'regionLayer',
               'MainfieldPlacement' : 'RegionLayer'})

lb.export_gpf(lb.merge([holes,regionLayer]), 
              {'FileName' : fn_prefix+'_holes.gpf',
               'GridResolution' : fine_res,
               'SubfieldSizeX' : fine_subfield_size,
               'SubfieldSizeY' : fine_subfield_size,
               'MainfieldSizeX' : fine_mainfield_size,
               'MainfieldSizeY' : fine_mainfield_size,
               'BeamStepSize' : fine_bss,
               'MainFieldResolution' : fine_mainfield_res,
               'SubFieldResolution' : fine_subfield_res,
               'FractureMode' : 'Curved',
               'RegionLayerSet' : 'regionLayer',
               'MainfieldPlacement' : 'RegionLayer'})

lb.export_gpf(lb.merge([phc,regionLayer]), 
              {'FileName' : fn_prefix+'_phc.gpf',
               'GridResolution' : fine_res,
               'SubfieldSizeX' : fine_subfield_size,
               'SubfieldSizeY' : fine_subfield_size,
               'MainfieldSizeX' : fine_mainfield_size,
               'MainfieldSizeY' : fine_mainfield_size,
               'BeamStepSize' : fine_bss,
               'MainFieldResolution' : fine_mainfield_res,
               'SubFieldResolution' : fine_subfield_res,
               # 'FractureMode' : 'Conventional',
               'FractureMode' : 'Curved',
               'RegionLayerSet' : 'regionLayer',
               'MainfieldPlacement' : 'RegionLayer'})
lb.export_gpf(negY,
              {'FileName' : fn_prefix+'_yJunction.gpf',
               'GridResolution' : fine_res,
               'SubfieldSizeX' : fine_subfield_size,
               'SubfieldSizeY' : fine_subfield_size,
               'MainfieldSizeX' : fine_mainfield_size,
               'MainfieldSizeY' : fine_mainfield_size,
               'BeamStepSize' : fine_bss,
               'MainFieldResolution' : fine_mainfield_res,
               'SubFieldResolution' : fine_subfield_res,
               'FractureMode' : 'LRFT',
               })

lb.export_gpf(phccover,
              {'FileName' : fn_prefix+'_regionLayer.gpf',
               'GridResolution' : fine_res,
               'SubfieldSizeX' : fine_subfield_size,
               'SubfieldSizeY' : fine_subfield_size,
               'MainfieldSizeX' : fine_mainfield_size,
               'MainfieldSizeY' : fine_mainfield_size,
               'BeamStepSize' : fine_bss,
               'MainFieldResolution' : fine_mainfield_res,
               'SubFieldResolution' : fine_subfield_res,
               'FractureMode' : 'Conventional',
               'RegionLayerSet' : 'regionLayer'})

### REIMPORT EXPORTED GPF FILES ###
print "Re-importing..."
holes = lb.import_gpf({'FileName' : fn_prefix+'_holes.gpf'})
phc = lb.import_gpf({'FileName' : fn_prefix+'_phc.gpf'})
yJunction = lb.import_gpf({'FileName' : fn_prefix+'_yJunction.gpf'})
tethers = lb.import_gpf({'FileName' : fn_prefix+'_tethers.gpf'})
fine = lb.import_gpf({'FileName' : fn_prefix+'_fine.gpf'})
coarse = lb.import_gpf({'FileName' : fn_prefix+'_coarse.gpf'})
# There is a trick here: the bias is used to make the regionlayer a
# rectangle again. Actually, to ensure that the crystal has no cuts, I will bias the rectangle to be larger than the crystal
regionLayer = lb.import_gpf({'FileName' : fn_prefix+'_regionLayer.gpf'})
regionLayer = lb.bias(regionLayer,{'Bias' : 30}) # to remove subfield cuts
regionLayer = lb.heal(regionLayer,{'TargetLayer' : 'regionLayer'}) # to ensure naming is good
fine = lb.merge([holes,phc,fine,tethers,yJunction])

### PEC ###
print "Performing PEC..."
pec_coarse = lb.heal(lb.grid(coarse,{'DatabaseGrid' : fine_res}),{'TargetLayer' : 'coarse'})
pec_merge = lb.merge([fine,pec_coarse])
pec = lb.pec(pec_merge,{'MaxNumOfDoseClasses' : 256,
                        'PSFFileName' : 'Y:\AndrewM\BEAMER\BEAMERpy\pec.xrz',
                        'PSFType' : 'Numeric'})
pec_fine = lb.extract(pec,{'LayerSet' : '0'})
pec_coarse = lb.extract(pec,{'LayerSet' : 'coarse'})

### EXPORTING AGAIN ###
print "Final export..."
lb.export_gpf(lb.merge([fine,regionLayer]),
              {'FileName' : fn_prefix+'_'+rotate_bln+'_fine.gpf',
               'GridResolution' : fine_res,
               'SubfieldSizeX' : fine_subfield_size,
               'SubfieldSizeY' : fine_subfield_size,
               'MainfieldSizeX' : fine_mainfield_size,
               'MainfieldSizeY' : fine_mainfield_size,
               'BeamStepSize' : fine_bss,
               'MainFieldResolution' : fine_mainfield_res,
               'SubFieldResolution' : fine_subfield_res,
               'FractureMode' : 'Conventional',
               'RegionLayerSet' : 'regionLayer',
               'MainfieldPlacement' : 'RegionLayer'})

lb.export_gpf(lb.merge([pec_fine,regionLayer]),
              {'FileName' : fn_prefix+'_'+rotate_bln+'_pec_fine.gpf',
               'GridResolution' : fine_res,
               'SubfieldSizeX' : fine_subfield_size,
               'SubfieldSizeY' : fine_subfield_size,
               'MainfieldSizeX' : fine_mainfield_size,
               'MainfieldSizeY' : fine_mainfield_size,
               'BeamStepSize' : fine_bss,
               'MainFieldResolution' : fine_mainfield_res,
               'SubFieldResolution' : fine_subfield_res,
               'FractureMode' : 'Conventional',
               'RegionLayerSet' : 'regionLayer',
               'MainfieldPlacement' : 'RegionLayer'})

lb.export_gpf(coarse,
              {'FileName' : fn_prefix+'_'+rotate_bln+'_coarse.gpf',
               'GridResolution' : coarse_res,
               'SubfieldSizeX' : coarse_subfield_size,
               'SubfieldSizeY' : coarse_subfield_size,
               'MainfieldSizeX' : coarse_mainfield_size,
               'MainfieldSizeY' : coarse_mainfield_size,
               'BeamStepSize' : coarse_bss,
               'MainFieldResolution' : coarse_mainfield_res,
               'SubFieldResolution' : coarse_subfield_res,
               'FractureMode' : 'Conventional'})

lb.export_gpf(pec_coarse,
              {'FileName' : fn_prefix+'_'+rotate_bln+'_pec_coarse.gpf',
               'GridResolution' : coarse_res,
               'SubfieldSizeX' : coarse_subfield_size,
               'SubfieldSizeY' : coarse_subfield_size,
               'MainfieldSizeX' : coarse_mainfield_size,
               'MainfieldSizeY' : coarse_mainfield_size,
               'BeamStepSize' : coarse_bss,
               'MainFieldResolution' : coarse_mainfield_res,
               'SubFieldResolution' : coarse_subfield_res,
               'FractureMode' : 'Conventional'})

### REMOVE INTERMEDIATE FILES ###
os.remove(fn_prefix+'_coarse.gpf')
os.remove(fn_prefix+'_fine.gpf')
os.remove(fn_prefix+'_phc.gpf')
os.remove(fn_prefix+'_holes.gpf')
os.remove(fn_prefix+'_tethers.gpf')
os.remove(fn_prefix+'_regionLayer.gpf')
os.remove(fn_prefix+'_yJunction.gpf')
