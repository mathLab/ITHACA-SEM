<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D MMF implicit diffusion </description>
    <executable>MMFDiffusion</executable>
    <parameters> MMFDiffPlane.xml</parameters>
    <files>
        <file description="Session File"> MMFDiffPlane.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-5">0.000110375</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5"> 0.00137303</value>
        </metric>
    </metrics>
</test>
