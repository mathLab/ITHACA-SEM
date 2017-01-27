<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Periodic 2D aerofoil mesh</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list 2d-periodic.mcf 2d-periodic.xml:xml:test</parameters>
    <files>
        <file description="Input File">2d-periodic.mcf</file>
        <file description="Input File 2">2d-periodic.STEP</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
