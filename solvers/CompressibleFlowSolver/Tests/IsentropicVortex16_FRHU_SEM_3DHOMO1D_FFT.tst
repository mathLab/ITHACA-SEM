<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=4, FRHU, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_FRHU_SEM_3DHOMO1D_FFT.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_FRHU_SEM_3DHOMO1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.00276e-06</value>
            <value variable="rhou" tolerance="1e-12">6.50207e-06</value>
            <value variable="rhov" tolerance="1e-12">6.10116e-06</value>
            <value variable="rhow" tolerance="1e-12">0</value>
            <value variable="E" tolerance="1e-12">2.20195e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">8.1427e-06</value>
            <value variable="rhou" tolerance="1e-12">2.15656e-05</value>
            <value variable="rhov" tolerance="1e-12">1.86306e-05</value>
            <value variable="rhow" tolerance="1e-12">0</value>
            <value variable="E" tolerance="1e-12">7.35728e-05</value>
        </metric>
    </metrics>
</test>


