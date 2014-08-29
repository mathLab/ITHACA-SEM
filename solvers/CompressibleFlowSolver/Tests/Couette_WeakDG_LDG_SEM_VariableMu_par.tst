<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, Variable Viscosity, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_SEM_VariableMu_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_SEM_VariableMu_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0809477</value>
            <value variable="rhou" tolerance="1e-12">51.9087</value>
            <value variable="rhov" tolerance="1e-8">0.221151</value>
            <value variable="E" tolerance="1e-12">4414.9</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0724519</value>
            <value variable="rhou" tolerance="1e-12">55.3223</value>
            <value variable="rhov" tolerance="1e-8">0.425157</value>
            <value variable="E" tolerance="1e-12">4226.24</value>
        </metric>
    </metrics>
</test>


