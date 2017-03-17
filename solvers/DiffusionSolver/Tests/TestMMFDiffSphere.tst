<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D MMF implicit diffusion </description>
    <executable>MMFDiffusion</executable>
    <parameters> TestMMFDiffSphere.xml</parameters>
    <files>
        <file description="Session File"> TestMMFDiffSphere.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-5">0.113019</value>
            <value variable="v" tolerance="1e-5">0.0437867</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">0.380462</value>
            <value variable="v" tolerance="1e-5">0.0905989</value>
        </metric>
    </metrics>
</test>
