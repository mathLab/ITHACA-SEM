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
            <value variable="rho" tolerance="1e-12">0.00406259</value>
            <value variable="rhou" tolerance="1e-12">0.00806001</value>
            <value variable="rhov" tolerance="1e-12">0.00754235</value>
            <value variable="E" tolerance="1e-12">0.0224819</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00462404</value>
            <value variable="rhou" tolerance="1e-12">0.0105546</value>
            <value variable="rhov" tolerance="1e-12">0.0107409</value>
            <value variable="E" tolerance="1e-12">0.0346782</value>
        </metric>
    </metrics>
</test>


