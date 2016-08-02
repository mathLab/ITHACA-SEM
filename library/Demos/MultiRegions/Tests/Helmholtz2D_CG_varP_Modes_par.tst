<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 2D CG with P=7, parallel</description>
    <executable>Helmholtz2D</executable>
    <parameters>Helmholtz2D_varP.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Helmholtz2D_varP.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-6">0.0045269</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.00420409</value>
        </metric>
    </metrics>
</test>


