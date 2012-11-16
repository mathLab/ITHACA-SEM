<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D CG, prisms, Neumann BCs, iterative ML, Par(3)</description>
    <executable>Helmholtz3D</executable>
    <parameters>helmholtz3D_Prism_Iter_ML_CG.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Helmholtz3D_Prism_Iter_ML_CG.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">2.24119e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">8.17542e-05</value>
        </metric>
    </metrics>
</test>
