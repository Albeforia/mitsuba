## Rendering Thin Transparent Layers with Extended Normal Distribution Functions
This is an implementation for the paper "Rendering Thin Transparent Layers with Extended Normal Distribution Functions".  
The code is based on Mitsuba renderer (http://www.mitsuba-renderer.org/) and
the BRDF plugin called *layeredconductor* for the new material is provided.

### Building
You can build the project using SCons. For details, please refer to the documentation of Mitsuba.  
Note if you want to build with CMake, the corresponding CMakeLists should be modified.

### Usage
Below is an example object with the material.

```xml
<shape type="sphere">
  <bsdf type="layeredconductor">
    <!-- concentration parameter of the upper surface -->
    <float name="kappa1" value="10000"/>

    <!-- concentration parameter of the lower surface -->
    <float name="kappa2" value="10"/>

    <!-- IOR of the transparent layer -->
    <float name="layerEta" value="1.5"/>
  </bsdf>
</shape>
```
