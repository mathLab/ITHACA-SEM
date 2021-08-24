<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D Helmholtz with boundary files, different Lagrange order</description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz3D_Cube_BndFile_Lagrange.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Cube_BndFile_Lagrange.xml</file>
        <file description="Session File">Helmholtz3D_Cube_BndFile_Lagrange.bc</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.59332e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.85419e-05</value>
        </metric>
    </metrics>
</test>
