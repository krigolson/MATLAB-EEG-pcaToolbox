ERP PCA Toolkit


Copyright (C) 1999-2018  Joseph Dien

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  Programs are not for medical use or any other application where use or misuse could cause injury.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

The rotations other than Varimax, Promax, and Infomax were based on code adapted from Matlab code made publicly available by Coen A. Bernaards and Robert I. Jennrich with their permission.
Website: http://www.stat.ucla.edu/research

The files ave_hdr_offsets_v.m, ave_hdr_offsets_v.m, ses_hdr_offsets_v.m, rd_onetr_allch.m, get_cell_offsets.m, get_fid, and wt_ses_hdr_v2.m are from EGI's EGI Toolbox with permission.  The files rd_PCAegis_hdr_v.m and wt_PCAave_hdr_v.m were modified from this same source.

The files el_topoplot.m, el_traditionaldipfit.m, and el_runica.m in the EEGlab folder in the External folder are renamed topoplot.m, traditionaldipfit.m, and runica.m files from EEGlab.  Only the names of the functions (and the commenting out of one line in runica) were changed and is provided under the GPL license.  This paragraph constitutes the required declaration that it has been modified from the original.  The topoplot modification was made so that calls could be made to the function without confusion occurring as to whether EEGlab's topoplot was intended or FieldTrip's topoplot was intended.  The traditionaldipfit modification was made so that the function would be available for use since it is not visible to outside programs under EEGlab’s standard installation procedure.  The commented out line in runica was performed so that output from repeated applications to a given dataset will fully replicate.

The file fmrib_fastr.m in the fmrib1.21 directory in the externals folder has two changes made to it to fix crashing errors.  The changes are documented in the code and included with the EP Toolkit, as permitted under GPL License.

The file amri_eeg_gac.m in the amri_eegfmri_toolbox directory in the externals folder has one set of changes made to it to fix crashing errors.  The changes are documented in the code and included with the EP Toolkit, as permitted under GPL License.

The file ft_sphericalSplineInterpolate.m in the fieldtrip folder in the External folder is unchanged from the fieldtrip-20140807 distribution and is included, as permitted under GPL license, as it is otherwise not directly accessible to EP Toolkit functions.

The file ft_sphericalSplineInterpolate.m in the fieldtrip folder in the External folder is a renamed sphericalSplineInterpolate.m file from the fieldtrip-20140807 distribution.  Only the names of the function was changed and is provided under the GPL license.  This paragraph constitutes the required declaration that it has been modified from the original.  The modification was made so that calls could be made to the function without confusion with the version accompanying fieldtrip.  This was done so that the function would be available for use since it is not visible to outside programs under fieldtrip’s standard installation procedure.

The files in the CCA directory in the externals folder are being included by permission of Maartens De Vos.  The compute_cca.m file has been modified to allow for options parameters to be passed on by the function call.

The files in the invcwt_v1 directory in the externals folder are being redistributed per permissions described in its license.txt file.