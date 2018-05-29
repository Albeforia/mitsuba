## A Physically-based Appearance Model for Special Effect Pigments
This is an implementation for our paper "A Physically-based Appearance Model for Special Effect Pigments".  
The code is based on Mitsuba renderer (http://www.mitsuba-renderer.org/) and the BRDF plugin called *glittery* for our model is provided.

### Building
You can build the project using SCons. For details, please refer to the documentation of Mitsuba.  
Note if you want to build with CMake, the corresponding CMakeLists should be modified.

### Usage
Below is an example object with our material.

```xml
<shape type="sphere">
    <bsdf type="glittery">
        <!-- roughness of the microfacet model -->
        <float name="alpha" value="0.04" />

        <!-- the complex IOR of the base layer -->
        <spectrum name="eta"  value="1.1" />
        <spectrum name="k"    value="1.5" />

        <!-- the IOR of the film -->
        <spectrum name="filmEta"  value="1.2" />

        <!-- thickness(nanometers) distribution of the film -->
        <spectrum name="height"      value="500" />
        <spectrum name="heightRange" value="50" />

        <!-- total number of facets -->
        <integer name="totalFacets" value="1000000"/>

        <float name="queryRadius" value="6"/>
    </bsdf>
</shape>
```
