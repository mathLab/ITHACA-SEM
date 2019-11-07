<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NavierStokes, restart from file, BCs from file, WeakDG, GLL</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters> hump3D_GLL.xml </parameters>
    <files>
        <file description="Session File">hump3D_GLL.xml </file>
        <file description="Session File">hump3D_b1.bc</file>
        <file description="Session File">hump3D_b2.bc</file>
        <file description="Restart File">hump3D.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">1.95793e-05</value>
            <value variable="rhou" tolerance="1e-12">0.0161709</value>
            <value variable="rhov" tolerance="1e-12">1.02091e-07</value>
             <value variable="rhow" tolerance="1e-12">0.000238819</value>
            <value variable="E" tolerance="1e-12">45.8424</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0608055</value>
            <value variable="rhou" tolerance="1e-12">55.621</value>
            <value variable="rhov" tolerance="1e-7">0.0616446</value>
            <value variable="rhow" tolerance="1e-12">0.444359</value>
            <value variable="E" tolerance="1e-12">60454.3</value>
        </metric>
    </metrics>
</test>