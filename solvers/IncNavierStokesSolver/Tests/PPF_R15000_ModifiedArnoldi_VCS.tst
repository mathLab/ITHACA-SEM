<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with Velocity Correction Scheme and a half mode  Ev = (2.4868e-03,1.5835e-01i) which corresponds to a multiplier of (1.00249,0.158351). Note without the restart file it need 1000's iterations</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters> PPF_R15000_3D.xml</parameters>
    <files>
        <file description="Session File">PPF_R15000_3D_VCS.xml</file>
        <file description="Session File">PPF_R15000_3D_VCS.rst</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value index="0" tolerance="0.001">1.00249,0.158351</value>
            <value index="1" tolerance="0.001">1.00249,-0.158351</value>
        </metric>
    </metrics>
</test>


