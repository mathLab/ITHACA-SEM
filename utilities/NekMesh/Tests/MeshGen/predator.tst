<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Complex geometry of predator drone, no bl</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list predator.mcf predator.xml:xml:test</parameters>
    <files>
        <file description="Input File">predator.mcf</file>
        <file description="Input File 2">predator.STEP</file>
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
