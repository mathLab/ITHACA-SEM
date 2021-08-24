<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>2D Helmholtz/Steady Diffusion With Direct Full Solve P=7</description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz2D_DirectFull.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_DirectFull.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.58061e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.80329e-05</value>
        </metric>
    </metrics>
</test>
