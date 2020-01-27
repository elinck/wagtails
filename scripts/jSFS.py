import allel
import moments
import numpy as np

# load vcf and extract folded sfs; write to file
vcf_file = "/data3/vcf_files/allopatric_parental/intermediate_files/demographic_inf_parentals_rename.vcf"
f = allel.read_vcf(vcf_file, fields='*')
gt = f['calldata/GT']
gt = allel.GenotypeArray(gt)
personata=gt[:,0:10]
alba=gt[:,10:20]
gn_p = personata.to_n_alt(fill=-1)
gn_a = alba.to_n_alt(fill=-1)
p_unlinked = allel.locate_unlinked(gn_p)
a_unlinked = allel.locate_unlinked(gn_a)
pass_linkage = np.logical_and(p_unlinked, a_unlinked)
final = gt.compress(pass_linkage, axis=0)
personata2=final[:,0:10]
alba2=final[:,10:20]
pers_ac = personata2.count_alleles()
alba_ac = alba2.count_alleles()
fsfs = allel.joint_sfs_folded(pers_ac, alba_ac)
m_fsfs = moments.Spectrum(fsfs)
m_fsfs.to_file("allopatric_parental_unlinked_sfs.txt")
