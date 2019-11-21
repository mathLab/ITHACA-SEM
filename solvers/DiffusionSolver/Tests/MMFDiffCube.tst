<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D MMF implicit diffusion </description>
    <executable>MMFDiffusion</executable>
    <parameters> MMFDiffCube.xml</parameters>
    <files>
        <file description="Session File"> MMFDiffCube.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-5">0.000339924</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5"> 0.00052942</value>
        </metric>
    </metrics>
</test>
