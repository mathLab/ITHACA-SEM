<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D MMF implicit diffusion </description>
    <executable>MMFDiffusion</executable>
    <parameters> MMFDiffSphere.xml</parameters>
    <files>
        <file description="Session File"> MMFDiffSphere.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-5">0.0142248</value>
            <value variable="v" tolerance="1e-5">0.00581111</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">0.0569858</value>
            <value variable="v" tolerance="1e-5">0.0270625</value>
        </metric>
    </metrics>
</test>
