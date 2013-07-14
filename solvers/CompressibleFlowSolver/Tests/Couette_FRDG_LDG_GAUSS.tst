<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRDG advection and LDG diffusion, GAUSS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRDG_LDG_GAUSS.xml</parameters>
    <files>
        <file description="Session File">Couette_FRDG_LDG_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0886025</value>
            <value variable="rhou" tolerance="1e-12">60.3537</value>
            <value variable="rhov" tolerance="1e-8">0.212262</value>
            <value variable="E" tolerance="1e-12">4866.89</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0738626</value>
            <value variable="rhou" tolerance="1e-12">61.0435</value>
            <value variable="rhov" tolerance="2e-6">0.22695</value>
            <value variable="E" tolerance="1e-12">4381.93</value>
        </metric>
    </metrics>
</test>


