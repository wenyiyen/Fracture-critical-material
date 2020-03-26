<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<h1 align="center">
    <a>Development of fracture critical material</a>
</h1>

<h4 align="center">
    <a href="#material-overview">Material Overview</a> |
     <a href="#material-compiling">Material Compiling</a> |
    <a href="#connection-model">Connection model</a> |
    <a href="#steel-frame-performance">Steel Frame Performance</a>
</h4>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->

This repo contains the descriptions of the fracture critical material and how to compile it onto Stanford Cluster System (Sherlock) for parallel computing. The use of the material is further described in the [connection model](#connection-model) and [steel frame performance](#steel-frame-performance) sections.

## Material Overview
In the proposed method, the connection is modeled using fiber beam elements, in which the hinge is discretized into fracture critical flanges and web fibers. In this model, the flange fiber constitutive (i.e. stress-strain) relationship represents the series of the weld and beam-hinging region of the connection. The flange fiber loses tensile resistance when fracture occurs, while the behavior under compressive forces is unaffected considering the re-contacting of the beam flange and column face. After connection fracture, the compressive resistance acts like a gap material, resembling the open and closing of the gap between the beam flange and column face.

* Material Constitutive<!--
<p align="center"><img style="max-width:500px" width="640" src="https://github.com/wenyiyen/Fracture-critical-material/blob/master/fig/constitutive.png" alt="constitutive"></p>-->
