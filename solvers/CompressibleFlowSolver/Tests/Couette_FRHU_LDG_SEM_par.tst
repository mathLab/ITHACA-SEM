<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRHU advection and LDG diffusion, SEM, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-metis Couette_FRHU_LDG_SEM_par.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Couette_FRHU_LDG_SEM_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0889265</value>
            <value variable="rhou" tolerance="1e-12">62.1008</value>
            <value variable="rhov" tolerance="1e-8">0.17592</value>
            <value variable="E" tolerance="1e-12">4904.01</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.076032</value>
            <value variable="rhou" tolerance="1e-12">56.1036</value>
            <value variable="rhov" tolerance="2e-6">0.265578</value>
            <value variable="E" tolerance="1e-12">4361.42</value>
        </metric>
    </metrics>
</test>


