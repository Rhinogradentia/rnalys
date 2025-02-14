Build image:

```
singularity build rnalys.sif rnalys.def

# or

apptainer build rnalys.sif rnalys.def
```

Run image:

```
apptainer run --bind 8050:8050 rnalys.sif

# or

singularity run --bind 8050:8050 rnalys.sif
```

Access App:

```
http://localhost:8050
```


Be aware! This does not use the original requirements.txt and install_packages.R --> this could be improved by importing them or using them directly of course.
