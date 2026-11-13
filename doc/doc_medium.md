# Medium Properties

The medium sets the sound speed, density and attenuation the simulated wave travels through, and the thermal properties that govern heating. Whether you simulate free water or a head, and which tissue properties you assume, therefore shapes the predicted focus, intensity and temperature rise. This page shows how to choose the medium and which tissue properties PRESTUS assumes by default.

On this page:

- [Choosing the Medium](#choosing-the-medium): the `water`, `layered` and `phantom` media and what each needs.
- [Default Tissue Properties](#default-tissue-properties): the acoustic and thermal defaults per tissue, and how much of the attenuation heats the tissue.

!!! note "At a glance"
    - **When you need it:** To choose between a free-water and a head simulation, or to change tissue properties (e.g. for a more conservative estimate).
    - **Before you start:** Nothing; the defaults apply until you change them.

## Choosing the Medium

Set `simulation.medium` to one of three values:

| Medium | What is simulated | Needs |
|---|---|---|
| `water` | The transducer in homogeneous water, placed at the inner edge of the requested PML layer. Free-water specifications and measurements refer to this. | No images |
| `layered` | A head built from the [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) segmentation. With the default `layers` and a [charm](https://simnibs.github.io/simnibs/build/html/documentation/command_line/charm.html) segmentation: water (including CSF), brain, skin, cortical skull bone and trabecular skull bone; the default `skull` layer has no charm label and is dropped. See [Layer mapping and smoothing](doc_preproc.md#layer-mapping-and-smoothing). The skull can instead take voxel-wise properties from a [(pseudo-)CT](doc_pseudoCT.md). | T1 (and T2/UTE), SimNIBS |
| `phantom` | A designed test object. The segmentation file is mapped onto medium masks without rotation, cropping or other [preprocessing](doc_preproc.md). | A `final_tissues.nii.gz` in a (dummy) SimNIBS output folder |

```yaml
simulation:
  medium: layered
```

Which layers a `layered` simulation contains is set under `layers`. A layer whose segmentation label is missing is dropped from the run. The HTML report lists the layers that were modelled and their properties.

A `phantom` [NIfTI](https://nifti.nimh.nih.gov/) holds artificial layers of a SimNIBS segmentation. It covers the simulation domain only: the PML is added outside the grid. Its voxels must be isotropic. The `phantom` medium accepts the same `layers` flexibility as `layered`, but applies no image processing, so 2D phantoms are a quick way to benchmark or debug.

## Default Tissue Properties

![medium_properties](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/medium_properties.png)

The properties below apply in a full layered simulation. They are set under `medium_properties` in `config_default.yaml` ([Parameters](doc_parameters.md#tissue-acoustic-thermal-properties-medium_properties)) and can be changed per tissue in your study config, for example for a more conservative estimate of intensity or heating. The defaults follow the literature and commercially available tools; they are not ground truth. The [uncertainty mode](doc_uncertainty.md) runs the [ITRUSST](https://itrusst.org/) bounds for you.

| Tissue | Density [kg/m³] | Sound speed [m/s] | Attenuation coefficient [dB/(cm·MHz^y)] | Attenuation power law (y) | Thermal conductivity [W/(m·K)] | Specific heat [J/(kg·K)] | Perfusion [mL/min/kg] | Absorption fraction | Equivalent absorption coefficient [dB/(cm·MHz^y)] |
|---|---|---|---|---|---|---|---|---|---|
| Water | 994 | 1500 | 0.00217 | 2 | 0.6 | 4178 | 0 | 1 | 0.00217 |
| Brain | 1046 | 1546 | 0.59 | 1.2 | 0.51 | 3630 | 559 | 1 | 0.59 |
| Skin | 1090 | 1610 | 0.4 | 1 | 0.37 | 3391 | 106 | 1 | 0.4 |
| Skull (single layer) | 1850 | 2800 | 13.3 | 1 | 0.32 | 1793 | 20 | 0.28 | 3.78 |
| Skull (cortical) | 1850 | 2800 | 13.3 | 1 | 0.32 | 1313 | 10 | 0.28 | 3.78 |
| Skull (trabecular) | 1700 | 2300 | 13.3 | 1 | 0.32 | 2274 | 30 | 0.28 | 3.78 |

The defaults set `absorption_fraction`; `absorption_coeff` stays empty. The last column is the coefficient to set instead (for skull, the [Pinton et al.](https://doi.org/10.1118/1.3668316) value from which the fraction 0.28 is rounded).

```yaml
# example: a more conservative skull for heating
medium_properties:
  skull_cortical:
    absorption_fraction: 0.5
```

### Absorption Fraction for Heating

For thermal simulations, PRESTUS uses only part of the attenuation as absorbed energy in the bioheat equation. Each tissue sets either `absorption_fraction` (a fraction of `alpha_coeff`) or `absorption_coeff` (a coefficient in dB/(cm·MHz^y) used directly); setting one empties the other ([Heat deposition](doc_simulations-thermal.md#heat-deposition)).

The skull default of 0.28 lies between the [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) default (16 %) and the [k-Plan](https://k-plan.io/) default (100 %). More measurements are needed to narrow the uncertainty about how much acoustic energy the skull absorbs.

??? note "Where the skull absorption fraction comes from"
    The bulk attenuation coefficient includes contributions from both longitudinal and shear waves. Longitudinal waves contribute most (about 80 %, [Wang et al., 2018](https://doi.org/10.1088/1674-1056/27/2/024302)) to absorption and heating in the skull (Pinton et al., 2012; [White et al., 2006](https://doi.org/10.1016/j.ultrasmedbio.2006.03.015)). For longitudinal and shear waves, absorption values of 2.7 dB/(cm·MHz) and 5.4 dB/(cm·MHz) have been reported (Pinton et al., 2012), with a bulk attenuation of 13.3 dB/(cm·MHz). PRESTUS therefore assumes absorption fraction = (longitudinal absorption + 0.2 × shear absorption) / bulk attenuation = (2.7 + 0.2 × 5.4) / 13.3 = 0.28. The equivalent `absorption_coeff` is 2.7 + 0.2 × 5.4 = 3.78 dB/(cm·MHz).

??? note "References"
    - Pinton, G., Aubry, J., Bossy, E., Muller, M., Pernot, M., & Tanter, M. Attenuation, scattering, and absorption of ultrasound in the skull bone. *Med Phys* **39**, 299–307 (2012). https://doi.org/10.1118/1.3668316
    - Wang, X.-D., Lin, W.-J., Su, C., & Wang, X.-M. Influence of mode conversions in the skull on transcranial focused ultrasound and temperature fields utilizing the wave field separation method: A numerical study. *Chin Phys B* **27**, 024302 (2018). https://doi.org/10.1088/1674-1056/27/2/024302
    - White, P. J., Clement, G. T., & Hynynen, K. Longitudinal and shear mode ultrasound propagation in human skull bone. *Ultrasound Med Biol* **32**, 1085–1096 (2006). https://doi.org/10.1016/j.ultrasmedbio.2006.03.015
