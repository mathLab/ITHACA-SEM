<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NavierStokes, restart from file, BCs from file, WeakDG, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters> hump3D_SEM.xml </parameters>
    <files>
        <file description="Session File">hump3D_SEM.xml </file>
        <file description="Session File">hump3D_b1.bc</file>
        <file description="Session File">hump3D_b2.bc</file>
        <file description="Restart File">hump3D.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">1.95847e-05</value>
            <value variable="rhou" tolerance="1e-12">0.0161757</value>
            <value variable="rhov" tolerance="1e-12">1.45248e-07</value>
             <value variable="rhow" tolerance="1e-12">0.000238819</value>
            <value variable="E" tolerance="1e-12">45.8422</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0602796</value>
            <value variable="rhou" tolerance="1e-12">55.5259</value>
            <value variable="rhov" tolerance="1e-12">0.0666846</value>
            <value variable="rhow" tolerance="1e-12">0.444353</value>
            <value variable="E" tolerance="1e-12">60454.4</value>
        </metric>
    </metrics>
</test>