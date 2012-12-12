<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady DG explicit diffusion, order 4, P=8</description>
    <executable>ADRSolver</executable>
    <parameters>ExDiffusion_m8.xml</parameters>
    <files>
        <file description="Session File">ExDiffusion_m8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.000337147</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00279079</value>
        </metric>
    </metrics>
</test>



