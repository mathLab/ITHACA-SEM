<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D Helmholtz/Steady Diffusion Reaction variable P with all element types and HDF5 output</description>
    <executable>ADRSolver</executable>
    <parameters>CubeAllElements_VarP.xml</parameters>
    <files>
        <file description="Session File">CubeAllElements_VarP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00130472</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00797808</value>
        </metric>
    </metrics>
</test>
