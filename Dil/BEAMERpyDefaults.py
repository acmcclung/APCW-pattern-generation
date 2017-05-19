
## all dimensions in microns ###
## IMPORT ##
#import_res =  0.0025
import_res =  0.0001
import_snap = 0.0001 
## EXPORT ##
## fine ##
fine_res = 0.0025 
fine_bss = 0.0025
fine_mainfield_res = 0.0025
fine_subfield_res = 0.00025
fine_mainfield_size = 160
fine_subfield_size = 4.075
## coarse ##
coarse_res = 0.01
coarse_bss = 0.01
coarse_mainfield_res = 0.01
coarse_subfield_res = 0.00047619
coarse_mainfield_size = 560
coarse_subfield_size = 4.525

## DICTIONARIES ##
extract_def = {'ExtentMode' : 'Maintain',
               'ExtractMode' : 'EntireLayoutExtract',
               'CellName' : '',
               'RegionLayer' : '',
               'RegionBehavior' : 'Clip',
               'AllExceptForRegions' : False,
               'DoseClassificationShrink' : 0,
               'ExtractBoxes' : []}
bias_def = {'SoftFrame' : 0.3, # microns of overlap between virtual tiles
            'CornerExtension' : 1, # described in 4.2.17.2
            'TargetLayer' : '0(0)',
            'Mode' : 'X-Y',
            'HierarchicalProcessing' : True, 
            # saves time if hierarchical, but may destroy hierarchy
            'LayerAssignment' : 'AllLayer',
            'SelectedLayerSet' : '*',
            'BiasLayerList' : []}
fracture_def = {# 'Resolution' : import_res,
    'BeamStepSize' : 1,
    'CurveTolerance' : 1,
    'FractureAxis' : 'X_AND_Y',
    'BssFracturing' : False,
    #                    'FractureAngle' : 'Rectangular',
    'FractureTolerance' : 1, # for 'curved' fracturing
    'FractureType' : 'Hierarchical',
    'FieldSizeX' : 1000, # for 'flat with fields'
    'FieldSizeY' : 1000,
    'SubfieldSizeX' : 0,
    'SubfieldSizeY' : 0,
    'TraversalDirection' : 'BottomUp',
    'FieldOverlapX' : 0, # MAYBE CHANGE THESE?
    'FieldOverlapY' : 0,
    'OverlapMethod' : 'Standard',
    'InterleavingSize' : 0,
    'InterlockLayer' : '*', # ???
    'MultipassMode' : 1,
    'MainfieldOffsetX' : 0,
    'MainfieldOffsetY' : 0,
    'SubfieldOffsetX' : 0,
    'SubfieldOffsetY' : 0,
    'MultipassLayer' : '*'}
fda_def = {'AssignmentType' : 'Assign'}
transform_def = {'CoordinateOrigin' : 'LayoutOrigin',
                 'ScaleX' : 1,
                 'ScaleY' : 1,
                 'ShiftX' : 0,
                 'ShiftY' : 0,
                 'ReflectX' : False,
                 'ReflectY' : False}
import_gpf_def = { 'MergeSubfieldCuts' : False }
export_gpf_def = { 'ExtentMode' : 'Maintain',
                   'LowerLeftX' : 0,
                   'LowerLeftY' : 0,
                   'UpperRightX' : 0,
                   'UpperRightY' : 0,
                   'FormatType' : '5000+ HS 100kV',
                   # 'GridResolution' : fine_res,
                   # 'MainFieldResolution' : fine_mainfield_res,
                   # 'SubFieldResolution' : fine_subfield_res,
                   # 'BeamStepSize' : fine_bss,
                   # 'MainfieldSizeX' : fine_mainfield_size,
                   # 'MainfieldSizeY' : fine_mainfield_size,
                   # 'SubfieldSizeX' : fine_subfield_size,
                   # 'SubfieldSizeY' : fine_subfield_size,
                   'CompactionRegionSize' : 10000, # what is this?
                   'DoseOrderingType' : 'AscendingDose',
                   'RegionTraversalMode' : 'MeanderX',
                   'FeatureOrderingType' : 'ArrayCompaction',
                   'SortedOrderLayer' : '*',
                   'MainfieldPlacement' : 'Fixed',
                   'YTrapezoids' : True,
                   'NormalizeToOne' : False,
                   'TrapezoidDoseCorrection' : False,
                   'ParallelogramCompaction' : True,
                   'SequenceFile' : '',
                   'FractureMode' : 'LRFT',
                   'RegionLayerSet' : 'regionLayer',
                   'BeamStepSizeFracturing' : True,
                   'SequenceLayer' : '',
                   'PlacementLayer' : '',
                   'NumberOfSleeves' : 1,
                   # 'AreaSelection' : 'RemainderWithSelected',
                   'AreaSelection' : 'SelectedThenFloating',
                   'FieldOverlapX' : 0,
                   'FieldOverlapY' : 0,
                   'OverlapMethod' : 'Standard',
                   'InterleavingSize' : 0,
                   'InterlockLayer' : '*',
                   'MultipassMode' : 1,
                   'MainfieldOffsetX' : 0,
                   'MainfieldOffsetY' : 0,
                   'SubfieldOffsetX' : 0,
                   'SubfieldOffsetY' : 0,
                   'MultipassLayer' : '*',
                   'UserDosePass' : [],
                   'OverlapMethod' : 'Standard',
                   'RegionList' : [] }
