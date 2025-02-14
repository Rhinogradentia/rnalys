Build container

```
docker build -t rnalys .
```

Run container

```
docker run -p 8050:8050 rnalys
```

Access app

```
http://localhost:8050
```

Be aware! This does not use the original requirements.txt and install_packages.R --> this could be improved by importing them or using them directly of course.
