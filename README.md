# auto_acpc_reorient : Cross-platform automatic AC-PC realignment/reorientation and coregistration robust to brain damage using SPM12

[![GitHub release](https://img.shields.io/github/release/lrq3000/auto_acpc_reorient.svg)](https://github.com/lrq3000/auto_acpc_reorient/releases/)

Cross-platform automatic AC-PC realignment/reorientation and coregistration for both healthy volunteers and brain damaged patients using template matching in SPM 12.

This is a set of routines to perform automatic reorient and coregistration with the toolbox [Statistical Parametric Mapping 12 (SPM12)](https://www.fil.ion.ucl.ac.uk/spm/).

![Automatic coregistration example using auto_acpc_coreg.m](img/coreg.png)

## Description

Setting up the AC-PC and reorienting images is a recurrent issue in between-subjects group analyses, since they rely on coregistration methods that, like the "unified segmentation" of SPM12, are for most sensitive to initial conditions (the starting orientation of the image).

The main function, `auto_acpc_reorient.m`, automatically (but approximately) calculates a reorientation transform onto a target template in MNI space, in two steps:

1. a non-linear coregistration of the input image onto a target template in MNI space is calculated using `spm_affreg`,
2. then another transform is calculated using Mutual Information on a joint histogram (spm_coreg), and then applies only the rigid-body transform part of both coregistrations to reorient the input image. This allows to mainly set the origin on the AC and correct for head rotation, in order to further proceed with the segmentation/normalisation of the image.

This whole reorientation scheme relies on the "template matching" principle (as in the SPM12 Old Normalize function), you therefore need to specify the appropriate template/reference image (we provide one, `t1group`, by default).

In any case, it is advised to visually check the automatically reoriented images afterwards, and [fix the orientation manually, using SPM -> Display](https://en.wikibooks.org/wiki/SPM/How-to#How_to_manually_change_the_orientation_of_an_image) if necessary (this process can be semi-automated with the [reorientation_registration_helper](https://github.com/lrq3000/pathmatcher#auxiliary-tool-reorientation-and-registration-helper) Python module).

Another function, `auto_acpc_coreg.m`, expands on the same ideas to allow coregistration between modalities (eg, functional BOLD on structural MPRAGE). It is advised that `auto_acpc_reorient()` to be first applied on the structural before applying `auto_acpc_coreg()` on the other modality (even if you do manually fix the reorientation, as this ensures that the T1 is somewhat in the MNI space, making it easier for `auto_acpc_coreg()` to find the correct translation matrix).

## Install

To install this tool :

* Add `spm12` to the path in MATLAB.
    * Optional: you can add `spm12\toolbox\OldNorm` to the path in MATLAB too if you want to use old affine method `affreg`, but we strongly disadvise it as it may cause reflections, even in rigid-body mode.
* Simply add this (auto_acpc_reorient) folder to the path in MATLAB too (or `cd` inside before typing the commands).

This will allow the commands `auto_acpc_reorient()` and `auto_acpc_coreg()` to be called from command-line (if no argument is given, a file selector dialog will open).

Note that this tool was tested on SPM12 v7771 (on 2024-06-27) and MATLAB R2023b.

## Usage

Very basic usage on a T1 (structural MRI) image:

```autoreorient('inputpath', 't1_orig.nii', 'regmode', 'jointhistogram')```

TODO: facade functions auto_acpc_reorient and auto_acpc_coreg to easily reorient or coregister images, using autoreorient in the background.

Helper facade functions are provided for convenience: `help auto_acpc_reorient`, for all the details and various options for reorientation of a T1. Type `help auto_acpc_coreg` for coregistering options.

Both facade scripts allows to use SPM filedialogs GUI, by simply typing `auto_acpc_reorient` or `auto_acpc_coreg` in the MATLAB prompt.

Note: by default, the 't1group' template will be used, which will use `T1_template_CAT12_rm_withskull.nii`. This is a template generated on 10 subjects using CAT12 that were manually reoriented to AC-PC and averaged, this provides better performance for reorientation than the more blurry MNI template. Note that this template is slightly better aligned to the AC-PC plane than the original MNI template, so that there may be a slight rotation bias compared to MNI if you use this custom template (usually it's mostly unnoticeable and this should have no influence if afterwards you do the SPM normalization on MNI on your data).

General guideline:

* If you want to reorient isotropic T1, use `auto_acpc_reorient`.
* If you want to reorient another modality (usually with less resolution), or an anisotropic T1, or an isotropic T1 but on your own custom T1 template, use `auto_acpc_coreg`.

Note that the scripts cannot be used from SPM12 GUI nor the BATCH system.

## Performance

There is no guarantee that this will work 100% of the times, although it was empirically observed to produce good results with our own data (young and old healthy or brain lesioned subjects of various severity up to extremely severe and of various etiologies including traumatic and anoxic, even with significant movement or metal artifacts or simulated nulling of whole brain areas).

The best results we got were by doing the following steps:

1. Call `auto_acpc_reorient()` on the structural in MNI space, so that it is also matching other SPM templates
2. Manually review and fix the misoriented structural images
3. Coregister the functional to SPM12 EPI template (this allows a correct translation matrix and a good basis for rotation)
4. Coregister the functional onto the structural (this fine-tunes rotation to precisely match the subject's structural)

The last 2 steps can be done by calling `auto_acpc_coreg()`, which has optimized default parameters for this task.

For a comparison of various methods for AC-PC reorientation, the following article is a very good read:

`Liu, Yuan, and Benoit M. Dawant. "Automatic detection of the anterior and posterior commissures on MRI scans using regression forests." 2014 36th Annual International Conference of the IEEE Engineering in Medicine and Biology Society. IEEE, 2014.`

## Citation

Please cite this work as following:

> auto_acpc_reorient, RRID:SCR_018308. [https://github.com/lrq3000/auto_acpc_reorient](https://github.com/lrq3000/auto_acpc_reorient)

## Authors

This software was developed by Stephen Karl Larroque (Coma Science Group, GIGA-Consciousness, University Hospital of Liege, Belgium).

The code is a rewrite from scratch, with several novel methods. Previous versions of auto_acpc_reorient.m were forks of [auto_reorient.m](https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;d1f675f1.0810) written by John Ashburner (FIL, UCL, London, UK) and Carlton Chu (FIL, UCL, London, UK), but not anymore, to allow to clarify/clean up the licensing (before, it was a [strange mix of Creative Commons with a discretionary possibility to relicense under GPLv3](https://en.wikibooks.org/wiki/SPM/How-to#How_to_automatically_reorient_images)).

## Frequently Asked Questions (FAQ)

* Can this script be used to "pseudo-normalize" the brain without using the usual unified segmentation workflow?

Theoretically, this should be possible by using a non-linear reorientation method, which is provided in SPM and in this script. A non-linear transform allows, on top of affine transforms, to apply independent resizing transforms in each of the 3 axes, hence allowing to resize the individual's brain to better match the template's and, supposedly, better match the brain areas.

However, non-linear transforms have one major issue: they can flip/mirror (ie, geometrical reflection) the brain, when the value of scaling is negative. This happens only when an odd number of axes scaling factors are negative (ie, 1 or 3 scaling factors will lead to a flipping, whereas 2 will not). There is currently no failsafe mechanism in SPM to prevent this from happening (and it's actually a quite complex issue mathematically), but more problematically, SPM does not detect this. The new experimental version of this project will include an autodetection feature to warn the user when a flipping/mirroring happened, and will even include some failsafe mechanisms to try to prevent this issue (but at the expense of approximating the non-linear reorientation result).

## Similar projects

We later found that K. Nemoto made a similar enhanced algorithm back in 2017, based on centering to the origin (instead of translating to match a template for us - what we call a pre-coregistration step) before applying the original auto_reorient.m script but tweaking it to apply a coregistration on a DARTEL template (instead of SPM12 ones in the original script, or a custom made template using CAT12 for this one). The resulting script, [acpc_coreg.m](https://web.archive.org/web/20180727093129/http://www.nemotos.net/scripts/acpc_coreg.m), can be found on [nemotos.net](https://www.nemotos.net/?p=1892).

Another earlier fork of auto_reorient.m named [spm_auto_reorient](https://github.com/CyclotronResearchCentre/spm_auto_reorient) was done by Christophe Phillips (Cyclotron Research Centre, University of Liege, Belgium), which itself forks the original by John Ashburner and Carlton Chu by adding support for more templates, option to apply the reorientation matrix to other files, sanitization of input variables, integration in the SPM12 main GUI using cfg files and usage in the BATCH system. An earlier version of `auto_acpc_reorient`, named `spm_auto_reorient_coregister`, was forked from `spm_auto_reorient`, before a major rewrite to rebase on the original codebase.

## License

The code is licensed under MIT Public License and authored by Stephen Karl Larroque except where indicated otherwise.

The templates are licensed under either MIT Public License or later, or CC0 (Public Domain), or the Unlicense at your convenience.

## Contact

If you have any feedback (questions, issues, suggestions), please contact Stephen Karl Larroque at: stephen dot larroque at uliege dot be.
