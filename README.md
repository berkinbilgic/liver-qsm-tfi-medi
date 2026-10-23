# Liver QSM: TFI and MEDI

MATLAB code and example in vivo data for abdominal (liver) QSM from 6-echo GRE (TE 1.02–8.12 ms, 2.89 T) with pilot-tone motion-resolved iMoCo reconstruction.

`script_abdominal_qsm_ref10_lbv_tfi_imoco_liverbbox_PT_cemre_v0.m` removes the background field with LBV and then runs two dipole inversions, both referenced to fat:

- **TFI-preconditioned LSQR**: no Laplacian in the forward model, edge-weighted gradient regularization and a fat-uniformity term (Velikina et al. 2023, Eq. 8, with the fat mean taken from χ itself), with the binary TFI preconditioner of Liu et al. 2017.
- **MEDI with SMV and fat regularization**: MEDI+0 with the fat mask in place of the CSF mask, as in Kang et al. 2023 (5 mm SMV, λ = 500, λ_fat = 100).

The last cell compares the two maps in the liver.

## Data (Git LFS)

The data in `data/` are stored with Git LFS. Install [git-lfs](https://git-lfs.com) before cloning, or run `git lfs pull` after cloning:

```
git clone https://github.com/berkinbilgic/liver-qsm-tfi-medi.git
```

GitHub's **Download ZIP** contains only the small LFS pointer files, not the data. Use `git clone`, or download each file in `data/` from its GitHub page.

- `exp100926_mgre_PT_imoco_softgating_fw.mat`: complex echoes (`img_cplx`), water, fat, R2* and field map (`fieldHz`)
- `exp100926_mgre_fb_coilsensmask.mat`: coil sensitivities (`receive`), used only for the body mask
- `liver_mask.mat`: TotalSegmentator (total_mr) liver mask, the union over the echo-2, RSS and water images, with its bounding box. The script uses it to normalize the data weight and for the liver comparison.

## Requirements

- MATLAB (tested with R2024b) and the Image Processing Toolbox
- Linux x86-64: the LBV solver is a compiled MEX file (`functions/mexMGv6.mexa64`)

## Usage

Run the script from the repository folder in MATLAB. The QSM maps are written to `qsm_output/` as NIfTI. MEDI writes its input file and `results/` to `medi_v2/` in the current folder.

## Included code

- `functions/`: TFI forward/adjoint operators, gradient operator, display and NIfTI helpers, dipole kernel, and LBV (`LBV.m`, `mexMGv6.mexa64`) from the MEDI toolbox
- `external/MEDI_toolbox/functions/`: `MEDI_L1` and its helpers, from the MEDI toolbox (Cornell MRI Research Lab)
- `external/NIfTI_20140122/`: NIfTI toolbox by Jimmy Shen (license in the folder)

## References

- Liu Z, et al. Preconditioned total field inversion (TFI) method for quantitative susceptibility mapping. MRM 2017;78:303–315.
- Velikina JV, et al. Data adaptive regularization with reference tissue constraints for liver quantitative susceptibility mapping. MRM 2023;90:385–399.
- Kang MS, et al. Free-breathing high isotropic resolution quantitative susceptibility mapping (QSM) of liver using 3D multi-echo UTE cones acquisition and respiratory motion-resolved image reconstruction. MRM 2023;90:1844–1858.
- Zhou D, et al. Background field removal by solving the Laplacian boundary value problem. NMR Biomed 2014;27:312–319.
- Liu Z, et al. MEDI+0: morphology enabled dipole inversion with automatic uniform cerebrospinal fluid zero reference for quantitative susceptibility mapping. MRM 2018;79:2795–2803.

Contact: berkin AT nmr.mgh.harvard.edu
