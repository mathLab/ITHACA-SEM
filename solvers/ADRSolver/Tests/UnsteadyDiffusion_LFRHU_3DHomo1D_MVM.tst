<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Homogeneous 1D Diffusion LFRHU MVM</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyDiffusion_LFRHU_3DHomo1D_MVM.xml</parameters>
    <files>
        <file description="Session File">UnsteadyDiffusion_LFRHU_3DHomo1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">3.71762e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">1.05767e-07</value>
        </metric>
    </metrics>
</test>
