<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D Helmholtz with boundary files</description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz3D_Cube_BndFile.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Cube_BndFile.xml</file>
        <file description="Session File">Helmholtz3D_Cube_BndFile.bc</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.53762e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.84523e-05</value>
        </metric>
    </metrics>
</test>
