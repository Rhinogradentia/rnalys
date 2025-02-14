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
