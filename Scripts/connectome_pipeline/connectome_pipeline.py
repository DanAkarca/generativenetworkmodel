#! /usr/bin/env python
import optparse
import os
import re
import sys


# ======================================================================
sys.path.append('/home/jb07/joe_python/GitHub/Modularity/')
sys.path.append('/home/jb07/python_modules/')

def main():
    p = optparse.OptionParser()

    p.add_option('--base_directory', '-b')
    p.add_option('--subject_list', '-s')
    p.add_option('--template_directory', '-t')
    p.add_option('--out_directory', '-o')
    p.add_option('--parcellation_directory', '-p')
    p.add_option('--acquisition_parameters', '-a')
    p.add_option('--index_file', '-i')
    sys.path.append(os.path.realpath(__file__))

    options, arguments = p.parse_args()
    base_directory = options.base_directory
    out_directory = options.out_directory
    subject_list = options.subject_list
    subject_list = [subject for subject in subject_list.split(
        ',') if subject]
    template_directory = options.template_directory
    parcellation_directory = options.parcellation_directory
    acquisition_parameters = options.acquisition_parameters
    index_file = options.index_file
    subjects_dir = out_directory + '/connectome/FreeSurfer/'

    if not os.path.isdir(out_directory + '/connectome/'):
        os.mkdir(out_directory + '/connectome/')
        os.mkdir(subjects_dir)

    os.environ['SUBJECTS_DIR'] = subjects_dir

    def connectome(subject_list, base_directory, out_directory):

        # ==================================================================
        # Loading required packages
        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as util
        from nipype.interfaces.freesurfer import ApplyVolTransform
        from nipype.interfaces.freesurfer import BBRegister
        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.diffusion_toolkit as dtk
        from nipype.interfaces.utility import Merge
        import numpy as np
        from additional_interfaces import AtlasValues
        from additional_interfaces import AparcStats
        from additional_interfaces import CalcMatrix
        from additional_interfaces import FreeSurferValues
        from additional_interfaces import Tractography
        from additional_pipelines import DWIPreproc
        from additional_pipelines import SubjectSpaceParcellation
        from additional_pipelines import T1Preproc

        from nipype import SelectFiles
        import os

        # ==================================================================
        # Defining the nodes for the workflow

        # Getting the subject ID
        infosource = pe.Node(interface=util.IdentityInterface(
            fields=['subject_id']), name='infosource')
        infosource.iterables = ('subject_id', subject_list)

        # Getting the relevant diffusion-weighted data
        templates = dict(T1='{subject_id}/anat/{subject_id}_T1w.nii.gz',
                         dwi='{subject_id}/dwi/{subject_id}_dwi.nii.gz',
                         bvec='{subject_id}/dwi/{subject_id}_dwi.bvec',
                         bval='{subject_id}/dwi/{subject_id}_dwi.bval')

        selectfiles = pe.Node(SelectFiles(templates),
                              name='selectfiles')
        selectfiles.inputs.base_directory = os.path.abspath(base_directory)

        # ==============================================================
        # T1 processing
        t1_preproc = pe.Node(interface=T1Preproc(), name='t1_preproc')
        t1_preproc.inputs.out_directory = out_directory + '/connectome/'
        t1_preproc.inputs.template_directory = template_directory

        # DWI processing
        dwi_preproc = pe.Node(interface=DWIPreproc(), name='dwi_preproc')
        dwi_preproc.inputs.out_directory = out_directory + '/connectome/'
        dwi_preproc.inputs.acqparams = acquisition_parameters
        dwi_preproc.inputs.index_file = index_file
        dwi_preproc.inputs.out_directory = out_directory + '/connectome/'

        # Eroding the brain mask
        erode_mask = pe.Node(interface=fsl.maths.ErodeImage(), name='erode_mask')

        # Reconstruction and tractography
        tractography = pe.Node(interface=Tractography(), name='tractography')
        tractography.iterables = ('model', ['CSA', 'CSD'])

        # smoothing the tracts
        smooth = pe.Node(interface=dtk.SplineFilter(
            step_length=0.5), name='smooth')

        # Moving to subject space
        subject_parcellation = pe.Node(interface=SubjectSpaceParcellation(), name='subject_parcellation')
        subject_parcellation.inputs.source_subject = 'fsaverage'
        subject_parcellation.inputs.source_annot_file = 'aparc'
        subject_parcellation.inputs.out_directory = out_directory + '/connectome/'
        subject_parcellation.inputs.parcellation_directory = parcellation_directory

        # Co-registering T1 and dwi
        bbreg = pe.Node(interface=BBRegister(), name='bbreg')
        bbreg.inputs.init='fsl'
        bbreg.inputs.contrast_type='t2'

        applyreg = pe.Node(interface=ApplyVolTransform(), name='applyreg')
        applyreg.inputs.interp = 'nearest'
        applyreg.inputs.inverse = True

        # Merge outputs to pass on to CalcMatrix
        merge = pe.Node(interface=Merge(3), name='merge')

        # calcuating the connectome matrix
        calc_matrix = pe.MapNode(interface=CalcMatrix(), name='calc_matrix', iterfield=['scalar_file'])
        calc_matrix.iterables = ('threshold', np.arange(0,20,10))

        # Getting values of diffusion measures
        FA_values = pe.Node(interface=AtlasValues(), name='FA_values')
        RD_values = pe.Node(interface=AtlasValues(), name='RD_values')
        AD_values = pe.Node(interface=AtlasValues(), name='AD_values')
        MD_values = pe.Node(interface=AtlasValues(), name='MD_values')

        # Getting additional surface measures
        aparcstats = pe.Node(interface=AparcStats(), name='aparcstats')
        aparcstats.inputs.parcellation_name = 'aparc'

        freesurfer_values = pe.Node(interface=FreeSurferValues(), name='freesurfer_values')
        freesurfer_values.inputs.parcellation_name = 'aparc'

        # ==================================================================
        # Setting up the workflow
        connectome = pe.Workflow(name='connectome')

        # Reading in files
        connectome.connect(infosource, 'subject_id', selectfiles, 'subject_id')

        # DWI preprocessing
        connectome.connect(infosource, 'subject_id', dwi_preproc, 'subject_id')
        connectome.connect(selectfiles, 'dwi', dwi_preproc, 'dwi')
        connectome.connect(selectfiles, 'bval', dwi_preproc, 'bvals')
        connectome.connect(selectfiles, 'bvec', dwi_preproc, 'bvecs')

        # CSD model and streamline tracking
        connectome.connect(dwi_preproc, 'mask', erode_mask, 'in_file')

        connectome.connect(selectfiles, 'bvec', tractography, 'bvec')
        connectome.connect(selectfiles, 'bval', tractography, 'bval')
        connectome.connect(dwi_preproc, 'dwi', tractography, 'in_file')
        connectome.connect(dwi_preproc, 'FA', tractography, 'FA')
        connectome.connect(erode_mask, 'out_file', tractography, 'brain_mask')

        # Smoothing the trackfile
        connectome.connect(tractography, 'out_track', smooth, 'track_file')

        # Preprocessing the T1-weighted file
        connectome.connect(infosource, 'subject_id', t1_preproc, 'subject_id')
        connectome.connect(selectfiles, 'T1', t1_preproc, 'T1')
        connectome.connect(t1_preproc, 'wm', subject_parcellation, 'wm')
        connectome.connect(t1_preproc, 'subjects_dir', subject_parcellation, 'subjects_dir')
        connectome.connect(t1_preproc, 'subject_id', subject_parcellation, 'subject_id')

        # Getting the parcellation into diffusion space
        connectome.connect(t1_preproc, 'subject_id', bbreg, 'subject_id')
        connectome.connect(t1_preproc, 'subjects_dir', bbreg, 'subjects_dir')
        connectome.connect(dwi_preproc, 'b0', bbreg, 'source_file')

        connectome.connect(dwi_preproc, 'b0', applyreg, 'source_file')
        connectome.connect(bbreg, 'out_reg_file', applyreg, 'reg_file')
        connectome.connect(subject_parcellation, 'renum_expanded', applyreg, 'target_file')

        # Calculating the FA connectome
        connectome.connect(tractography, 'out_track', calc_matrix, 'track_file')
        connectome.connect(dwi_preproc, 'FA', merge, 'in1')
        connectome.connect(dwi_preproc, 'RD', merge, 'in2')
        connectome.connect(tractography, 'GFA', merge, 'in3')
        connectome.connect(merge, 'out', calc_matrix, 'scalar_file')
        connectome.connect(applyreg, 'transformed_file', calc_matrix, 'ROI_file')

        # Getting values for additional measures
        connectome.connect(dwi_preproc, 'FA', FA_values, 'morpho_filename')
        connectome.connect(dwi_preproc, 'RD', RD_values, 'morpho_filename')
        connectome.connect(dwi_preproc, 'AD', AD_values, 'morpho_filename')
        connectome.connect(dwi_preproc, 'MD', MD_values, 'morpho_filename')
        connectome.connect(applyreg, 'transformed_file', FA_values, 'atlas_filename')
        connectome.connect(applyreg, 'transformed_file', RD_values, 'atlas_filename')
        connectome.connect(applyreg, 'transformed_file', AD_values, 'atlas_filename')
        connectome.connect(applyreg, 'transformed_file', MD_values, 'atlas_filename')

        # Getting FreeSurfer morphological values
        connectome.connect(t1_preproc, 'subject_id', aparcstats, 'subject_id')
        connectome.connect(t1_preproc, 'subjects_dir', aparcstats, 'subjects_dir')
        connectome.connect(aparcstats, 'lh_stats', freesurfer_values, 'lh_filename')
        connectome.connect(aparcstats, 'rh_stats', freesurfer_values, 'rh_filename')

        # ==================================================================
        # Running the workflow
        connectome.base_dir = os.path.abspath(out_directory)
        connectome.write_graph()
        connectome.run()

    os.chdir(out_directory)
    connectome(subject_list, base_directory, out_directory)

if __name__ == '__main__':
    # main should return 0 for success, something else (usually 1) for error.
    sys.exit(main())
