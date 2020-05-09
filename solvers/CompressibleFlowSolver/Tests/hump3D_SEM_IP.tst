<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NavierStokes, restart from file, BCs from file, WeakDG, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters> hump3D_SEM_IP.xml </parameters>
    <files>
        <file description="Session File">hump3D_SEM_IP.xml </file>
        <file description="Session File">hump3D_b1.bc</file>
        <file description="Session File">hump3D_b2.bc</file>
        <file description="Restart File">hump3D.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">1.95847e-05</value>
            <value variable="rhou" tolerance="1e-12">1.61757e-02</value>
            <value variable="rhov" tolerance="1e-12">1.76848e-07</value>
             <value variable="rhow" tolerance="1e-12">2.38819e-04</value>
            <value variable="E" tolerance="1e-12">4.58422e+01</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">6.02796e-02</value>
            <value variable="rhou" tolerance="1e-12">5.55254e+01</value>
            <value variable="rhov" tolerance="1e-12">8.03285e-02</value>
            <value variable="rhow" tolerance="1e-12">4.44353e-01</value>
            <value variable="E" tolerance="1e-12">6.04544e+04</value>
        </metric>
    </metrics>
</test>