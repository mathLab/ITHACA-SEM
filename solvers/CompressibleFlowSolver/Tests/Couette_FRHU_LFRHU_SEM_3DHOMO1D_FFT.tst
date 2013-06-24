<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRHU advection and LFRHU diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRHU_LFRHU_SEM_3DHOMO1D_FFT.xml</parameters>
    <files>
        <file description="Session File">Couette_FRHU_LFRHU_SEM_3DHOMO1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000397736</value>
            <value variable="rhou" tolerance="1e-12">48.1307</value>
            <value variable="rhov" tolerance="1e-8">0.145848</value>
            <value variable="rhow" tolerance="1e-8">9.12914e-06</value>
            <value variable="E" tolerance="1e-12">17520.2</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00140074</value>
            <value variable="rhou" tolerance="1e-12">83.3598</value>
            <value variable="rhov" tolerance="1e-8">0.505645</value>
            <value variable="rhow" tolerance="1e-8">2.94422e-05</value>
            <value variable="E" tolerance="1e-12">18956.6</value>
        </metric>
    </metrics>
</test>


