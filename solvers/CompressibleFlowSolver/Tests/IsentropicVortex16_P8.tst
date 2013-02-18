<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=8</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P8.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">2.77895e-08</value>
            <value variable="rhou" tolerance="1e-12">7.39606e-08</value>
            <value variable="rhov" tolerance="1e-12">1.19731e-07</value>
            <value variable="E" tolerance="1e-12">4.11526e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.23085e-07</value>
            <value variable="rhou" tolerance="1e-12">2.24846e-07</value>
            <value variable="rhov" tolerance="1e-12">5.57644e-07</value>
            <value variable="E" tolerance="1e-12">1.32341e-06</value>
        </metric>
    </metrics>
</test>


