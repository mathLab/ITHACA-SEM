<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=3 GAUSS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P3_GAUSS.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P3_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.00474302</value>
            <value variable="rhou" tolerance="1e-12">0.00947318</value>
            <value variable="rhov" tolerance="1e-12">0.00869843</value>
            <value variable="E" tolerance="1e-12">0.0254323</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00539172</value>
            <value variable="rhou" tolerance="1e-12">0.0113971</value>
            <value variable="rhov" tolerance="1e-12">0.0117187</value>
            <value variable="E" tolerance="1e-12">0.0374887</value>
        </metric>
    </metrics>
</test>


